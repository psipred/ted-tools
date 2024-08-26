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
		cout << " Usage: " << arg << " dist.mat" <<" seq.dat"<< " cont.mat" << endl;
	exit(EXIT_SUCCESS);
}
int main(int argc,const char* argv[])
{
	if (argc < 4) {
		print_help(argv[0]);
	}
	//count row number
	fstream fin(argv[1], ios::in);//open file
	char c;
	int length = 0;
	while (fin.get(c)) {
		if (c == '\n') {
			length++;
		}
	}
	//Get Distance Matrix
	double** DistanceMatrix;
	NewArray(&DistanceMatrix, length, length);
	freopen(argv[1], "r", stdin);
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < length; j++) {
			cin>> DistanceMatrix[i][j];
		}
	}
	fin.close();
	// Get pro Matrix
	double** ProMatrix;
	NewArray(&ProMatrix, length, length);
	freopen(argv[3], "r", stdin);
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < length; j++) {
			cin>> ProMatrix[i][j];
		}
	}
	fin.close();

	double** ContactMatrix;
	NewArray(&ContactMatrix, length, length);
	Cal_Con_Matrix_seq(ContactMatrix, DistanceMatrix, length,ProMatrix);
	vector<int> vsec;
	/***********************seq.dat****************************/
	ifstream ss(argv[2]);
	string row;
	while (getline(ss, row)) {
		vsec.push_back(atof(row.substr(15, 1).c_str()));
	}
    ss.close();
	/*********************************************************************************/
	/*                               up-bottom: cut domains                          */
	/*********************************************************************************/
	vector<int> W;
	for (int i = 0; i < length; i++) {
		W.push_back(i);
	}
	vector<string> N(length);
	for (int i = 0; i <length; i++) {
		string n = to_string(i+1);
		N[i] = n;
	}
	vector<string> initialFrag;
	int t = 0;
	string domainStr = Cut_Domain_sec(ContactMatrix, W, N,vsec);
	if (string::size_type id = domainStr.find(";") != string::npos) {
		//find the subdomain
		vector<string> domain;
		Split_String(domainStr, domain, ";");
		for (int j = 0; j < domain.size(); j++) {
			initialFrag.push_back(domain[j]);
		}
	}
	else {
		initialFrag.push_back(domainStr);
	}
	do {
		t = 0;
		vector<string> finalFrag;
		for (unsigned int j = 0; j < initialFrag.size(); j++)
		{
			vector<int> W1,V1;
			vector<string> N1;
			vector<string> spiltFrag;
			W1 = string_to_vector_int(initialFrag[j], W, N);
			N1 = string_to_vector_str(initialFrag[j], N, N);
			V1 = string_to_vector_int(initialFrag[j], vsec, N);
			string domain_str = Cut_Domain_sec(ContactMatrix, W1, N1,V1);
			if (string::size_type id = domain_str.find(";") != string::npos) {
				Split_String(domain_str, spiltFrag, ";");
				for (unsigned int i = 0; i < spiltFrag.size(); i++) {
					finalFrag.push_back(spiltFrag[i]);
				}
				t += 1;
			}
			else {
				finalFrag.push_back(domain_str);
			}
		}
		initialFrag.swap(finalFrag);
	} while (t > 0);
	/*for (int i = 0; i < initialFrag.size(); i++) {
		cout << initialFrag[i] << endl;
	}*/

	///*********************************************************************************/
	///*                               bottom-up: merge domains                        */
	///*********************************************************************************/
	vector<string> combine_domain = Combine_Domain(initialFrag, W, N, ContactMatrix);
	/*********************************************************************************/
	/*                               bottom-up: fliter                               */
	/*********************************************************************************/
	vector<string> domainsOut;
	if (combine_domain.size() > 1) {
		int q = 0;
		do {
			for (int i = 0; i < combine_domain.size(); i++) {
				replace_All(combine_domain[i], ";", ",");
			}
			string combineFrag;
			vector<string> MergeDomain;
			for (int i = 0; i < combine_domain.size(); i++) {
				double max_nnc = 0.0;
				vector<int> W1 = string_to_vector_int(combine_domain[i], W, N);
				double W1nnc = Ave_Density(ContactMatrix, W1) / 2;
				//cout << combine_domain[i] << ":" << W1nnc << endl;
				if (W1nnc < 1) {
					for (int j = 0; j < i; j++) {
						vector<int> W2 = string_to_vector_int(combine_domain[j], W, N);
						//double W2nnc = Ave_Density(ContactMatrix, W2) / 2;
						double nnc_W1_W2 = Inter_nnc(ContactMatrix, W1, W2);
						if (max_nnc < nnc_W1_W2) {
							max_nnc = nnc_W1_W2;
							combineFrag = combine_domain[i] + ";" + combine_domain[j];
						}
					}
					for (int j = i + 1; j < combine_domain.size(); j++) {
						vector<int> W2 = string_to_vector_int(combine_domain[j], W, N);
						//double W2nnc = Ave_Density(ContactMatrix, W2) / 2;
						double nnc_W1_W2 = Inter_nnc(ContactMatrix, W1, W2);
						if (max_nnc < nnc_W1_W2) {
							max_nnc = nnc_W1_W2;
							combineFrag = combine_domain[i] + ";" + combine_domain[j];
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
				for (int i = 0; i < combine_domain.size(); i++) {
					if (is_string_in_vector(spiltSegment, combine_domain[i]) != 1) {
						MergeDomain.push_back(combine_domain[i]);
					}
				}
				combine_domain.swap(MergeDomain);
				q = combine_domain.size();
			}
			else {
				domainsOut.swap(combine_domain);
				q = 0;
			}
		} while (q > 0);
	}
	else {
		domainsOut.swap(combine_domain);
	}
	///*********************************************************************************/
	///*                               Output Domain                                   */
	///*********************************************************************************/
	vector<string> domain;
	for (int i = 0; i < domainsOut.size(); i++) {
		vector<string> subDomain;
		Split_String(domainsOut[i], subDomain, ",");
		if (subDomain.size() == 1) {
			vector<string> c;
			Split_String(domainsOut[i], c, "~");
			string s = trim(c[0]) + "~" + trim(c[1]);
			domain.push_back(s);
		}
		else {
			vector<int> Domain_N_num;
			vector<int> Miss_N_num;
			for (int i = 0; i < subDomain.size(); i++) {
				vector<int> N1_num;
				N1_num = string_to_vector_int(subDomain[i], W, N);
				Domain_N_num.insert(Domain_N_num.end(), N1_num.begin(), N1_num.end());
			}
			sort(Domain_N_num.begin(), Domain_N_num.end());
			string s;
			vector<int> Dis_continue;
			vector<string> frag;
			for (int t = 0; t < Domain_N_num.size() - 1; t++) {
				int pro = Domain_N_num[t], suf = Domain_N_num[t + 1];
				int profix = int(atof(N[pro].substr(0, 4).c_str()));
				int suffix = int(atof(N[suf].substr(0, 4).c_str()));
				if (suffix == profix + 1 || suffix == profix) {
					continue;
				}
				else if (suffix != profix + 1) {
					Dis_continue.push_back(t);
				}
			}
			if (Dis_continue.size() != 0) {
				string sf = trim(N[Domain_N_num[0]]) + "~" + trim(N[Domain_N_num[Dis_continue[0]]]);
				frag.push_back(sf);
				int i = 0;
				if (Dis_continue.size() > 1) {
					do {
						int start = Dis_continue[i] + 1, end = Dis_continue[i + 1];;
						string sm = trim(N[Domain_N_num[start]]) + "~" + trim(N[Domain_N_num[end]]);
						frag.push_back(sm);
						i++;
					} while (i < Dis_continue.size() - 1);
				}
				string se = trim(N[Domain_N_num[Dis_continue.back() + 1]]) + "~" + trim(N[Domain_N_num.back()]);
				frag.push_back(se);
				for (int i = 0; i < frag.size() - 1; i++) {
					s += (frag[i] + ",");
				}
				s += frag.back();
			}
			else {
				s = trim(N[Domain_N_num[0]]) + "~" + trim(N[Domain_N_num.back()]);
			}
			domain.push_back(s);
		}
	}
	string str;
	for (int i = 0; i < domain.size() - 1; i++) {
		str += (domain[i] + "/");
	}
	str += domain.back();
	cout << str << endl;
	return 0;

}