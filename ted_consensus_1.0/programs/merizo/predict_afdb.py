from __future__ import annotations

import os
import argparse
import textwrap
import zipfile
import time

from scipy.optimize import linear_sum_assignment

import torch

from model.network import Merizo

from model.utils.features import generate_features_domain
from model.utils.utils import (
    get_device,
    format_dom_str,
    instance_matrix,
    write_pdf_predictions,
    write_pdb_predictions,
    write_fasta,
    clean_domains,
    clean_singletons,
    get_ids,
    remap_ids,
    shuffle_ids,
    separate_components,
)

# --- Constants for cleanup and iterative segmentation

MIN_DOMAIN_SIZE = 25    # minimum number of residues in a domain
MIN_FRAGMENT_SIZE = 5  # minimum number of residues in a single segment
DOM_AVE = 250           # half of the average domain size of CATH / for iteration mode
CONF_THRESHOLD = 0.95   # minimum domain confidence / for iteration mode

script_dir = os.path.dirname(os.path.abspath(__file__))

def iterative_segmentation(
        network: torch.nn.Module, 
        inputs: tuple[torch.tensor, torch.tensor, torch.tensor, torch.tensor, torch.tensor], 
        domain_ids: torch.tensor, 
        conf_res: torch.tensor,
        max_iterations: int,
    ) -> tuple[torch.tensor, torch.tensor]:
    """_summary_

    Args:
        network (torch.nn.Module): _description_
        inputs (tuple[torch.tensor, torch.tensor, torch.tensor, torch.tensor, torch.tensor]): _description_
        domain_ids (torch.tensor): _description_
        conf_res (torch.tensor): _description_

    Returns:
        tuple[torch.tensor, torch.tensor]: _description_
    """
    
    n_iterations = 0
    iterate = True
    ignore_index = []
    
    while iterate:
        ids, _ = get_ids(domain_ids)
        unique_ids = {}
        
        for d in ids:
            if d.item() not in ignore_index:
                dom_mask = domain_ids == d
                dom_conf_res = conf_res[dom_mask]
                dom_nres = len(dom_conf_res)
                
                assert len(dom_conf_res.unique()) == 1
                
                dom_conf = dom_conf_res.mean()
                
                cond1 = dom_nres > DOM_AVE # * 2
                # cond2 = dom_conf < CONF_THRESHOLD
                
                if cond1: # and cond2:
                    unique_ids[d.item()] = dom_conf
                else:
                    ignore_index.append(d.item())
                    
        if len(unique_ids) == 0:
            iterate = False
            break
        
        counter = 1
        for k, _ in unique_ids.items():
            domain_mask = domain_ids == k
            
            domain_ids_, conf_res_ = network(inputs, mask=domain_mask)
            
            ids, ndoms_ = get_ids(domain_ids_)

            # Get new confidence scores per new domain
            new_dom_conf = torch.zeros_like(ids).float()
            for i, d in enumerate(ids):
                new_mask = domain_ids_ == d
                new_dom_conf[i] = conf_res_[new_mask].mean()
                
            # If the segment is still 1 domain, skip
            if ndoms_ == 1:
                ignore_index.append(k)
                
            # Otherwise split and overwrite old domain ids and confidences
            else:
                # Offset new ids by at least no_classes to ensure no overlap
                # with the old ids, and then transplant back in
                dd = domain_ids_ + torch.max(domain_ids)
                dd[domain_ids_ == 0] = 0

                domain_ids[domain_mask] = dd
                
                # Assign all residues in the new domain with the predicted
                # domain confidence
                conf_res[domain_mask] = conf_res_
                counter += 1
                
        n_iterations += 1
        if n_iterations == max_iterations:
            iterate = False
            break
        
    return domain_ids, conf_res, n_iterations

def read_split_weight_files(directory):
    weights = {}
    
    # Get the list of weight files in the directory
    weight_files = [file for file in os.listdir(directory) if file.endswith('.pt')]
    
    # Read the weights from each file
    for file in weight_files:
        file_path = os.path.join(directory, file)
        subset_weights = torch.load(file_path, map_location=lambda storage, loc: storage)
        
        # Add the subset weights to the overall weights dictionary
        weights.update(subset_weights)
    
    return weights

def merge_doms(domain_ids, dm, ri, max_merge, d0=8.0, d=1.5, alpha=0.43, beta=0.95, ri_sep=4, min_s=0.3, min_d=8.0):
    """ Calculate the domain interaction score for pairs of domains following 
        the UniDoc method.
        
        d0, d, alpha, beta, ri_sep are hyperparameters reported in Zhu et al., 
        Bioinformatics, 39(2), 2023, bta070.
    
    """
    
    domain_ids_ = domain_ids.clone()

    # Contact probability based on CA-CA distance
    # Zhu et al., Equation (1)
    p = 1 / (1 + torch.exp((dm - d0) / d))
    r = ri[:, :, None] - ri[:, None, :] > ri_sep
    pr = p * r.int()
    
    run_merge = True
    n_merge = 0
    
    while run_merge:
        dom_ids = torch.unique(domain_ids_[domain_ids_.nonzero()])
        S = torch.zeros(len(dom_ids), len(dom_ids))

        # Calculate the interaction score between residues of the same domain
        # Zhu et al., Equation (3)
        dis_intra = {}
        for d in dom_ids:
            l = 1 / (len(domain_ids_[domain_ids_ == d]) ** beta)
            d_idx = domain_ids_ == d
            dp = pr[:, d_idx][:, :, d_idx]
            d_intra = l * torch.sum(dp)
            dis_intra[d.item()] = d_intra

        dom_combs = torch.combinations(dom_ids, 2)
        
        # Calculate the interaction score between residues of pairs of domains
        # Zhu et al., Equation (2)
        for dij in dom_combs:
            di, dj = dij
            
            di_mask = domain_ids_ == di
            dj_mask = domain_ids_ == dj
            
            dij_dist = torch.min(dm[:, di_mask][:, :, dj_mask])
            
            if dij_dist <= min_d:
                d_intra = torch.min(dis_intra[di.item()], dis_intra[dj.item()])
                
                li = len(domain_ids_[di_mask])
                lj = len(domain_ids_[dj_mask])

                l = 1 / ((li ** alpha) * (lj ** alpha))
                d_inter = l * torch.sum(p[:, domain_ids_ == di][:, :, domain_ids_ == dj])
                sij = d_inter - d_intra
            else:
                sij = 0.

            si = (dom_ids == di).nonzero().item()
            sj = (dom_ids == dj).nonzero().item()
            
            S[si, sj] = S[sj, si] = sij # make symmetric otherwise LSA doesn't work in the next step

        S_ = S.clone()
        S_[S_ < min_s] = 0.

        if torch.any(S > 0):
            assigned = []
            
            # Get most optimal pairs of assignments:
            mask = S_.sum(dim=-1) != 0
            assign_ids = dom_ids.clone().detach()[mask]
            assignable = S_[mask]

            row_ind, col_ind = linear_sum_assignment(-assignable)
            
            # For each domain, find and merge with the most positive domain in
            new_dom_ids = domain_ids_.clone()
            for r, c in zip(row_ind, col_ind):
                if dom_ids[c] not in assigned and assign_ids[r] not in assigned:
                    new_dom_ids[new_dom_ids == dom_ids[c]] = assign_ids[r]
                    # print("Merged domain {} into domain {}".format(dom_ids[c], assign_ids[r]))
                    
                    assigned.extend([dom_ids[c], assign_ids[r]])

            # Overwrite existing domain ids
            domain_ids_ = new_dom_ids
            n_merge += 1
        else:
            run_merge = False
            
        if n_merge == max_merge:
            run_merge = False
    
    return domain_ids_, n_merge

def segment(network, args, pdb_path, device, zipped, outfile):
    
    start_time = time.time()
    
    if zipped:
        pdb_path = pdb_path.filename
        
    pdb_name = os.path.basename(pdb_path)
    outname = pdb_path[:-4] + "_merizo"
    
    if args.outdir is not None:
        outname = os.path.join(args.outdir, os.path.basename(outname))
                
    s, z, r, t, ri, pdb, _, md5 = generate_features_domain(pdb_path, device, zipped=zipped)
    nres = s.shape[1]

    inputs = (s, z, r, t, ri)
    domain_ids, conf_res = network(inputs)

    # Resegment without the background residues
    # -------------------------------------------
    
    nres_domain = domain_ids.count_nonzero()
    mask = domain_ids != 0.

    if args.dual_pass and torch.any(mask) == 1 and nres_domain > MIN_DOMAIN_SIZE:
        dr_ids, dr_conf = network(inputs, mask=mask)
        
        domain_ids[mask] = dr_ids
        conf_res[mask] = dr_conf
    
    # -------------------------------------------
    
    nres_domain = domain_ids.count_nonzero()
    nres_ndomain = nres - nres_domain

    # If --iterate mode, iteratively segment domains 
    n_iterations = 0
    
    if args.iterate:
        if nres_domain > DOM_AVE * 2:
            domain_ids, conf_res, n_iterations = iterative_segmentation(
                network, inputs, domain_ids, conf_res, args.max_iterations)

    R_pred = instance_matrix(domain_ids)[0]
    domain_ids = separate_components(R_pred, z, domain_ids)

    if len(torch.unique(domain_ids)) > 1:
        domain_ids = clean_domains(domain_ids, MIN_DOMAIN_SIZE)
        domain_ids = clean_singletons(domain_ids, MIN_FRAGMENT_SIZE)
        
    if args.merge:
        # Merge based on UniDoc algorith
        domain_ids, n_merge = merge_doms(domain_ids, z.squeeze(-1), ri, args.max_merge)
        
        R_pred = instance_matrix(domain_ids)[0]
        domain_ids = separate_components(R_pred, z, domain_ids)

    # Recompute the domain map given the new assignment
    R_pred = instance_matrix(domain_ids)[0]
    
    if args.shuffle_indices:
        domain_ids = shuffle_ids(domain_ids)
    else:
        domain_ids = remap_ids(domain_ids)
    
    conf_global = conf_res.mean()
    
    _, ndoms = get_ids(domain_ids)

    # --------------
    # Outputs 
    # --------------
    
    dom_str = format_dom_str(domain_ids, ri)
    
    # Calculate again
    nres = domain_ids.shape[0]
    # nres_domain = domain_ids.count_nonzero()
    # nres_ndomain = nres - nres_domain

    if args.save_pdb or args.save_domains:
        write_pdb_predictions(
            pdb,
            domain_ids,
            conf_res,
            ri,
            args.save_domains,
            args.conf_filter,
            args.plddt_filter,
            outname,
        )

    if args.save_fasta:
        write_fasta(pdb, outname, pdb_name[:-4])

    if args.save_pdf:
        R_pred = instance_matrix(domain_ids)[0]
        p_conf = torch.sqrt(conf_res[None, :] * conf_res[:, None])
        p_conf = p_conf * R_pred

        title = "{} | {} predicted domains".format(pdb_name, ndoms)
        write_pdf_predictions(R_pred, p_conf, title, outname)

    if dom_str == '':
        dom_str = 'NULL'
        
    end_time = time.time() - start_time 
    bn, _ = os.path.splitext(pdb_name)
    
    with open(outfile, 'a+') as fn:
        fn.write("{}\t{}\t{}\t{}\t{}\t{:.5f}\n".format(
            bn,
            md5,
            nres, 
            ndoms, 
            dom_str,
            conf_global.item(),
            # n_iterations,
            # end_time,
            # device_name,
        ))


def main():
    # Read the config file
    parser = argparse.ArgumentParser(
        prog='ProgramName',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
            If you use Merizo, please cite the following paper:
                Lau, et al., 2023. Merizo: a rapid and accurate domain segmentation method using invariant point attention. bioRxiv, doi: https://doi.org/10.1101/2023.02.19.529114
            
            Example usage:
                python predict.py -d cpu -i examples/2xdqA.pdb
                python predict.py -d cpu -i examples/*.pdb --save_domains --save_pdf --save_fasta
                python predict.py -d cpu -i examples/2xdqA.pdb --save_domains --plddt_filter
                
            For AlphaFold2 models, the iterative segmentation routine may give better results on longer models:
                python predict.py -d cpu -i examples/AF-Q96PD2-F1-model_v4.pdb --iterate --plddt_filter 60 --conf_filter 0.75
         ''')
    )
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--input", type=str, nargs="+", required=False, help="Pass a single file or list of filenames as -i 1ubq.pdb or *.pdb.")
    group.add_argument("-l", "--list", type=str, required=False, help="Pass a file containing paths to files.")
    group.add_argument("-z", "--zipfile", type=str, required=False, help="Pass a zip archive containing member pdb files.")

    parser.add_argument("-d", "--device", type=str, default="cuda", help="Hardware to run on. Options: 'cpu', 'cuda', 'mps'.")
    parser.add_argument("-o", "--outdir", type=str, default=None, required=False, help="Directory to save models to.")
    parser.add_argument("--out", dest="outfile", type=str, required=True, help="File to save results to.")
    parser.add_argument("--save_pdf", dest="save_pdf", action="store_true", help="Include to save the domain map as a pdf.")
    parser.add_argument("--save_pdb", dest="save_pdb", action="store_true", help="Include to save the result as a pdb file. All domains will be included unless --conf_filter or --plddt_filter is used.")
    parser.add_argument("--save_domains", dest="save_domains", action="store_true", help="Include to save parsed domains as separate pdb files. Also saves the full pdb.")
    parser.add_argument("--save_fasta", dest="save_fasta", action="store_true", help="Include to save a fasta file of the input pdb.")
    parser.add_argument("--conf_filter", dest="conf_filter", type=float, default=None, help="(float, [0-1]) If specified, only domains with a pIoU above this threshold will be saved. ")
    parser.add_argument("--plddt_filter", dest="plddt_filter", type=float, default=None, help="(float, [0-1]) If specified, only domain with a plDDT above this threshold will be saved. Note: if used on a non-AF structure, this will correspond to crystallographic b-factors.")
    parser.add_argument("--iterate", dest="iterate", action="store_true", default=True, help=f"If used, domains under a length threshold (default: {DOM_AVE} residues) will be re-segmented.")
    parser.add_argument("--merge", dest="merge", action="store_true", default=True, help=f"If used, domains will be merged if inter-domain contacts are improved. Uses UniDoc's merge protocol.")
    parser.add_argument("--max_iterations", dest="max_iterations", type=int, default=3, help="(int [1, inf]) Specify the maximum number of re-segmentations that can occur.")
    parser.add_argument("--max_merge", dest="max_merge", type=int, default=3, help="(int [1, inf]) Specify the maximum number of re-segmentations that can occur.")
    parser.add_argument("--dual_pass", dest="dual_pass", action="store_true", default=True, help=f"If used, segment a second time without NDRs predicted in the first pass.")
    parser.add_argument("--shuffle_indices", dest="shuffle_indices", action="store_true", help="Shuffle domain indices - increases contrast between domain colours in PyMOL.")
    parser.add_argument("--label", dest="label", required=False, default='merizo-v2', help="Label to add to generated files.")
    parser.add_argument("--inherit_chopping", dest="inherit_chopping", required=False, default=None, help="NOT IMPLEMENTED")
    
    args = parser.parse_args()

    # args.input = ['examples/human_200_test/AF-O95302-F1-model_v4.pdb', 'examples/human_200_test/AF-O94779-F1-model_v4.pdb'] #, 'examples/AF-Q8NGE8-F1-model_v4_c2m.pdb']
    # args.input = ['examples/human_200_test/AF-P33527-F1-model_v4.pdb']
    # args.save_pdb = True
    # args.merge = True
    # args.iterate = True
    # args.outdir = 'examples/human_200_test/'
    # args.max_iterations = 2

    flag_count = sum([args.input is not None, args.list is not None, args.zipfile is not None])
    if flag_count != 1:
        parser.error('Exactly one flag (-i, -l or -z) with a file path must be provided.')
        
    if args.save_domains and args.outdir == None:
        parser.error("Need to provide -o/--outdir if --save_pdb is passed.")
        
    if args.outdir is not None:
        if not os.path.exists(args.outdir):
            os.mkdir(args.outdir)
    
    if args.input is not None:
        files = args.input

    if args.list is not None:
        with open(args.list, 'r') as f:
            files = [line.rstrip('\n') for line in f]
            
    if args.zipfile is not None:
        zf = zipfile.ZipFile(args.zipfile)
        files = [f for f in zf.filelist if os.path.basename(f.filename)[-4:] == '.pdb']
        zipped = zf
    else:
        zipped = False
        
    device = get_device(args.device)
    network = Merizo().to(device)

    network.load_state_dict(read_split_weight_files(os.path.join(script_dir, 'weights')), strict=True)
    network.eval()

    failed = []
    for pdb_path in files:
        try:  
            with torch.no_grad():
                segment(network, args, pdb_path, device, zipped=zipped, outfile=args.outfile)
        except:
            print("Failed: ", pdb_path)
            # continue
            if args.device == 'cuda':
                failed.append(pdb_path)

    if args.device == 'cuda': # Re-try failed models on CPU
        for pdb_path in failed:
            device = torch.device("cpu")
            network = network.to(device)

            try: 
                with torch.no_grad():
                    segment(network, args, pdb_path, device, zipped=zipped, outfile=args.outfile)
            except:
                print("{}\tSegmentation failed even on CPU".format(os.path.basename(pdb_path)))
        
        

if __name__ == "__main__":
    main()
