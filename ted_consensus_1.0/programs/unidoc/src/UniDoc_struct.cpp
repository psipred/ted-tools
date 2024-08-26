#include"utils.h"
#include <stdlib.h>
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include<vector>
#include <typeinfo>
#include<iterator>
#include <algorithm>
using namespace std;
void print_help(const char* arg) {
	cout << endl
		<< " Usage: " << arg << " PDB.pdb" << "  chain" << endl << endl;
	exit(EXIT_SUCCESS);
}
int main(int argc, const char* argv[])
{
	if (argc < 2) {
		print_help(argv[0]);
	}
	//step 1:get PDB file
	string pdb = argv[1];
	ifstream inPDB(argv[1]);
	string PDBstr;
	int tag = 1;
	while (getline(inPDB, PDBstr)) {
		string chain = argv[2];
		if (PDBstr.compare(21, 1, chain) == 0) {
			
			if (PDBstr.substr(0, 5) == "ATOM ") {
				if (PDBstr.substr(12, 4) == " CA ") {
					tag = 1;
				}
				if (PDBstr.substr(12, 4) != " CA ") {
					tag = 0;//have secondary structure
					break;
				}
			}
			if (PDBstr.substr(0, 3) == "TER") break;
		}
	}
	inPDB.close();
	/*********************************************************************************/
	/*                               get atom number                                 */
	/*********************************************************************************/
	int length = 0;
	ifstream inPDB1(argv[1]);
	string PDBstr1;
	string residueNum1 = "-10001";
	while (getline(inPDB1, PDBstr1)) {
		string chain = argv[2];
		if (PDBstr1.compare(21, 1, chain) == 0) {
			if (PDBstr1.substr(0, 5) == "ATOM ") {
				string atom_name = PDBstr1.substr(12, 4);
				string atom_num = PDBstr1.substr(22, 5);
				if (atom_name == " CA " || atom_name == " CB ") {
					if (residueNum1 != atom_num) {
						residueNum1 = atom_num;
						length++;
					}
				}
			}
			if (PDBstr1.substr(0, 3) == "TER") break;
		}
	}
	inPDB1.close();
	/*********************************************************************************/
	/*                               get CB coordinate                               */
	/*********************************************************************************/
	double X0[length], Y0[length], Z0[length];
	int N0[length];
	int m = -1, flag = 0;
	ifstream inPDB2(argv[1]);
	string PDBstr2;
	string residueNum2 = "-1000";
	while (getline(inPDB2, PDBstr2))
	{
		string chain = argv[2];
		if (PDBstr2.compare(21, 1, chain) == 0) {
			if (PDBstr2.substr(0, 5) == "ATOM ") {
				string atom_name = PDBstr2.substr(12, 4);
				string amin_abbr = PDBstr2.substr(17, 3);
				string atom_num = PDBstr2.substr(22, 5);
				if (atom_name == " CA " || atom_name == " CB ") {
					if (residueNum2 != atom_num) {
						residueNum2 = atom_num;
						flag = 0;
						m++;
					}
				}
				if (atom_name == " CB " || (flag == 0 && atom_name == " CA ")) {
					X0[m] = atof(PDBstr2.substr(30, 8).c_str());
					Y0[m] = atof(PDBstr2.substr(38, 8).c_str());
					Z0[m] = atof(PDBstr2.substr(46, 8).c_str());
					N0[m] = int(atof(PDBstr2.substr(22, 4).c_str()));
					if (atom_name == " CB ")flag = 1;
				}
			}
			if (PDBstr2.substr(0, 3) == "TER") break;
		}
	}
	inPDB2.close();
	vector<double> X(X0, X0 + length);
	vector<double> Y(Y0, Y0 + length);
	vector<double> Z(Z0, Z0 + length);
	vector<int> N(N0, N0 + length);
	vector<int> W;
	for (int i = 0; i < N.size(); i++) {
		W.push_back(i);
	}
	/*********************************************************************************/
	/*                     Calculate the Contact Matrix                              */
	/*********************************************************************************/
	double** ContactMatrix;
	NewArray(&ContactMatrix, length, length);
	Cal_Con_Matrix(ContactMatrix, X, Y, Z);
	/*********************************************************************************/
    /*                               get secondary structure                         */
    /*********************************************************************************/
	vector<string> ss;
	vector<int> vsec;
	if (tag == 0) {
		string chain = argv[2];
		string s = "./stride " + std::string("-r") + chain + " " + pdb + ">pdb_ss";
		system(s.c_str());
		ifstream fin("pdb_ss");
		string row;
		while (getline(fin, row)) {
			if (row.substr(0, 3) == "ASG") {
				string text = row.substr(14, 1);
				char a[1];
				a[0] = text[0];
				if (isupper(int(a[0]))) {
					continue;
				}
				else {
					ss.push_back(row);
				}
			}
		}
		fin.close();
		int k = 0;
		for (int i = 0; i < ss.size(); i++) {
			if (ss[i].substr(0, 3) == "ASG") {
				string SecType = ss[i].substr(24, 1);
				if (SecType == "C" || SecType == "T") {
					vsec.push_back(1);
				}
				else if (SecType == "B") {
					vsec.push_back(2);
				}
				else {
					vsec.push_back(3);
				}
			}
		}
	}
	/*********************************************************************************/
	/*                               up-bottom: cut domains                          */
	/*********************************************************************************/
	vector<string> initialFrag;
	vector<string> domainsOut;
	int t = 0;
	if (tag == 0) {
		if (string::size_type id = Cut_Domain_sec(ContactMatrix, W, N, vsec).find(";") != string::npos) {
			//find the subdomain
			vector<string> domain;
			Split_String(Cut_Domain_sec(ContactMatrix, W, N, vsec), domain, ";");
			for (int j = 0; j < domain.size(); j++) {
				initialFrag.push_back(domain[j]);
			}
		}
		else {
			initialFrag.push_back(Cut_Domain_sec(ContactMatrix, W, N, vsec));
		}
		do {
			t = 0;
			vector<string> finalFrag;
			for (unsigned int j = 0; j < initialFrag.size(); j++)
			{
				vector<int> W1, N1, V1;
				vector<string> spiltFrag;
				W1 = string_to_vector(initialFrag[j], W, N);
				N1 = string_to_vector(initialFrag[j], N, N);
				V1 = string_to_vector(initialFrag[j], vsec, N);
				if (string::size_type id = Cut_Domain_sec(ContactMatrix, W1, N1, V1).find(";") != string::npos) {
					Split_String(Cut_Domain_sec(ContactMatrix, W1, N1, V1), spiltFrag, ";");
					for (unsigned int i = 0; i < spiltFrag.size(); i++) {
						finalFrag.push_back(spiltFrag[i]);
					}
					t += 1;
				}
				else {
					finalFrag.push_back(Cut_Domain_sec(ContactMatrix, W1, N1, V1));
				}
			}
			initialFrag.swap(finalFrag);
		} while (t > 0);
	}
	else {
		if (string::size_type id = Cut_Domain(ContactMatrix, W, N).find(";") != string::npos) {
			//find the subdomain
			vector<string> domain;
			Split_String(Cut_Domain(ContactMatrix, W, N), domain, ";");
			for (int j = 0; j < domain.size(); j++) {
				initialFrag.push_back(domain[j]);
			}
		}
		else {
			initialFrag.push_back(Cut_Domain(ContactMatrix, W, N));
		}
		do {
			t = 0;
			vector<string> finalFrag;
			for (unsigned int j = 0; j < initialFrag.size(); j++)
			{
				vector<int> W1, N1;
				vector<string> spiltFrag;
				W1 = string_to_vector(initialFrag[j], W, N);
				N1 = string_to_vector(initialFrag[j], N, N);
				if (string::size_type id = Cut_Domain(ContactMatrix, W1, N1).find(";") != string::npos) {
					Split_String(Cut_Domain(ContactMatrix, W1, N1), spiltFrag, ";");
					for (unsigned int i = 0; i < spiltFrag.size(); i++) {
						finalFrag.push_back(spiltFrag[i]);
					}
					t += 1;
				}
				else {
					finalFrag.push_back(Cut_Domain(ContactMatrix, W1, N1));
				}
			}
			initialFrag.swap(finalFrag);
		} while (t > 0);
	}
	/*********************************************************************************/
	/*                               bottom-up: merge domains                        */
	/*********************************************************************************/
	int q = 0;
	do {
		for (int i = 0; i < initialFrag.size(); i++) {
			replace_All(initialFrag[i], ";", ",");
		}
		q = 0;
		vector<string> MergeDomain;
		string combineFrag;
		double MaxNnc = 0.0;
		for (unsigned int i = 0; i < initialFrag.size(); i++) {
			//The first sub_domain
			vector<int> W1, N1;
			W1 = string_to_vector(initialFrag[i], W, N);
			N1 = string_to_vector(initialFrag[i], N, N);
			for (int j = 0; j < i; j++) {
				//The second sub_domain
				vector<int> W2, N2;
				W2 = string_to_vector(initialFrag[j], W, N);
				N2 = string_to_vector(initialFrag[j], W, N);
				//Step 3.1:calculate the inter_nnc and intra_nnc
				//THE Inter nnc
				double nnc_W1_W2 = Inter_nnc(ContactMatrix, W1, W2);
				//THE Intra nnc i.e the domain's average density
				double W1nnc = Ave_Density(ContactMatrix, W1) / 2;
				double W2nnc = Ave_Density(ContactMatrix, W2) / 2;
				double MINnnc = min(W1nnc, W2nnc);
				//Judge the two domain whether can be merged
				if (nnc_W1_W2 >= MINnnc) {
					double dif_nnc = nnc_W1_W2 - MINnnc; //Amount of change
					if (MaxNnc < dif_nnc) {
						MaxNnc = dif_nnc;
						combineFrag = initialFrag[i] + ";" + initialFrag[j];
					}
				}
			}
		}
		vector<string> segement, spiltSegment;
		if (combineFrag.size() != 0) {
			MergeDomain.push_back(combineFrag);
			Split_String(MergeDomain[0], segement, ";");
			string addDomain = segement.back();
			spiltSegment.push_back(addDomain);
			int a = addDomain.size();
			int b = MergeDomain[0].size();
			string beginDomain = MergeDomain[0].substr(0, b - a - 1);
			spiltSegment.push_back(beginDomain);
			for (int i = 0; i < initialFrag.size(); i++) {
				if (is_string_in_vector(spiltSegment, initialFrag[i]) != 1) {
					MergeDomain.push_back(initialFrag[i]);
				}
			}
			initialFrag.swap(MergeDomain);
			q = initialFrag.size();
		}
		else {
			domainsOut.swap(initialFrag);
			q = 0;
		}
	} while (q > 0);
	/*for (int i = 0; i < domainsOut.size(); i++) {
		cout << domainsOut[i] << endl;
	}*/
	vector<string> domain;
	for (int i = 0; i < domainsOut.size(); i++) {
		vector<string> subDomain;
		Split_String(domainsOut[i], subDomain, ",");
		if (subDomain.size() == 1) {
			domain.push_back(domainsOut[i]);
		}
		else {
			vector<int> Domain_N;
			vector<int> Miss_N;
			for (int i = 0; i < subDomain.size(); i++) {
				vector<int> N1;
				N1 = string_to_vector(subDomain[i], N, N);
				Domain_N.insert(Domain_N.end(), N1.begin(), N1.end());
			}
			sort(Domain_N.begin(), Domain_N.end());
			string s;
			for (int i = N[0]; i <= N[N.size() - 1]; i++) {
				if (is_element_in_vector(Domain_N, i) != 1) {
					Miss_N.push_back(i);
				}
			}
			/*for (int i = 0; i < Miss_N.size(); i++) {
				cout << Miss_N[i] << endl;
			}*/
			if (Miss_N.size() != 0) {
				int k = 0;
				vector<string> frag;
				if (Miss_N.size() > 1) {
					if (Domain_N[0] < Miss_N[0]) {
						string sf = to_string(Domain_N[0]) + "~" + to_string(Miss_N[0] - 1);
						frag.push_back(sf);
					}
					do {
						if (Miss_N[k + 1] != Miss_N[k] + 1) {
							string sm = to_string(Miss_N[k] + 1) + "~" + to_string(Miss_N[k + 1] - 1);
							frag.push_back(sm);
						}
						k++;
					} while (k < Miss_N.size() - 1);

					if (Domain_N.back() > Miss_N.back()) {
						string se = to_string(Miss_N.back() + 1) + "~" + to_string(Domain_N.back());
						frag.push_back(se);
					}
				}
				else {
					if (Domain_N[0] > Miss_N[0]) {
						string sf = to_string(Domain_N[0]) + "~" + to_string(Domain_N.back());
						frag.push_back(sf);
					}
					else if(Domain_N.back() < Miss_N[0]) {
						string sf = to_string(Domain_N[0]) + "~" + to_string(Domain_N.back());
						frag.push_back(sf);
					}
					else {
						string sm1 = to_string(Domain_N[0]) + "~" + to_string(Miss_N[0]-1);
						frag.push_back(sm1);
						string sm2 = to_string(Miss_N[0] + 1) + "~" + to_string(Domain_N.back());
						frag.push_back(sm2);
					}
				}
				for (int i = 0; i < frag.size() - 1; i++) {
					s += (frag[i] + ",");
				}
				s += frag.back();
			}
			else {
				s = to_string(Domain_N[0]) + "~" + to_string(Domain_N.back());
			}
			domain.push_back(s);
		}
	}
	string str;
	for (int i = 0; i < domain.size()-1; i++) {
		str += (domain[i] + "/");
	}
	str += domain.back();
	cout << str << endl;
	return 0;
}