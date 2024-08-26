#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>
#include <string.h>
#include <malloc.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cassert>
using namespace std;
int lc, uc, lb, ub;
/*
  delete " "
*/
string trim(string& s) {
	if (!s.empty()) {
		s.erase(0, s.find_first_not_of(" "));
		s.erase(s.find_last_not_of(" ") + 1);
	}
	return s;
}
/*
Judege the string whether in vector
para:
	vector<string> v
	int e
*/
bool is_string_in_vector(vector<string> v, string e) {
	vector<string>::iterator it;
	it = find(v.begin(), v.end(), e);
	if (it != v.end()) {
		return true;
	}
	else {
		return false;
	}
}
/*
Spilt the string according to string c
para:
	const string& s-- String need to be cut
	vector<string>& v--the vector<string> after cutting
	const string& c--According to the cut string
return:
	vector<string>& v
*/
void Split_String(const string& s, vector<string>& v, const string& c)
{
	string::size_type pos1, pos2;
	pos2 = s.find(c);
	pos1 = 0;
	while (string::npos != pos2) {
		v.push_back(s.substr(pos1, pos2 - pos1));
		pos1 = pos2 + c.size();
		pos2 = s.find(c, pos1);
	}
	if (pos1 != s.length()) {
		v.push_back(s.substr(pos1));
	}
}
/*
String replacement
para:
	string& s--the string
	string& old_value--the string need to be change
	string& new_value-- the string use to change
return:
	string s
*/
string& replace_All(string& s, const string& old_value, const string& new_value)
{
	for (string::size_type pos(0); pos != string::npos; pos += new_value.length()) {
		if ((pos = s.find(old_value)) != string::npos) {
			s.replace(pos, old_value.length(), new_value);
		}
		else { break; }
	}
	return s;
}
/*
Transform the string to vector
para:
	string s--the string need to be changed
	vector<string> W--the vector of initial domain
	vector<int> N--the Residue serial number of  initial domain
return:
	vector<string> W1 ---the vector of string
*/
vector<string> string_to_vector_str(const string s, const vector<string> W, const vector<string> N)
{
	vector<string> W1;
	vector<string> domainparser;
	if (string::size_type tag = s.find(";") != string::npos) {
		Split_String(s, domainparser, ";");
		for (unsigned int t = 0; t < domainparser.size(); t++) {
			vector<string> b, split_frag;
			if (string::size_type id = domainparser[t].find(",") != string::npos) {
				Split_String(domainparser[t], b, ",");
				for (unsigned int k = 0; k < b.size(); k++) {
					vector<string> c;
					Split_String(b[k], c, "~");
					string c0 = c[0], c1 = c[1];
					for (int n = 0; n < N.size(); n++) {
						if (c0 == N[n])  lc = n;
						if (c1 == N[n])  uc = n;
					}
					for (int i = lc; i <= uc; i++) {
						W1.push_back(W[i]);
					}
				}
			}
			else {
				Split_String(domainparser[t], b, "~");
				string b0, b1;
				b0 = b[0], b1 = b[1];
				for (int n = 0; n < N.size(); n++) {
					if (b0 == N[n])lb = n;
					if (b1 == N[n])ub = n;
				}
				for (int i = lb; i <= ub; i++) {
					W1.push_back(W[i]);
				}
			}
		}
	}
	else {
		vector<string> b, split_frag;
		if (string::size_type id = s.find(",") != string::npos) {
			Split_String(s, b, ",");
			for (unsigned int k = 0; k < b.size(); k++) {
				vector<string> c;
				Split_String(b[k], c, "~");
				string c0, c1;
				c0 = c[0], c1 = c[1];
				for (int n = 0; n < N.size(); n++) {
					if (c0 == N[n])lc = n;
					if (c1 == N[n])uc = n;
				}
				for (int i = lc; i <= uc; i++) {
					W1.push_back(W[i]);
				}
			}
		}
		else {
			Split_String(s, b, "~");
			string b0, b1;
			b0 = b[0], b1 = b[1];
			for (int n = 0; n < N.size(); n++) {
				if (b0 == N[n]) lb = n;
				if (b1 == N[n]) ub = n;
			}
			for (int i = lb; i <= ub; i++) {
				W1.push_back(W[i]);
			}
		}
	}
	return W1;
}
/*
Transform the string to vector
para:
	string s--the string need to be changed
	vector<int> W--the vector of initial domain
	vector<int> N--the Residue serial number of  initial domain
return:
	vector<int> W1 ---the vector of string
*/
vector<int> string_to_vector_int(const string s, const vector<int> W, const vector<string> N)
{
	vector<int> W1;
	vector<string> domainparser;
	if (string::size_type tag = s.find(";") != string::npos) {
		Split_String(s, domainparser, ";");
		for (unsigned int t = 0; t < domainparser.size(); t++) {
			vector<string> b, split_frag;
			if (string::size_type id = domainparser[t].find(",") != string::npos) {
				Split_String(domainparser[t], b, ",");
				for (unsigned int k = 0; k < b.size(); k++) {
					vector<string> c;
					Split_String(b[k], c, "~");
					string c0 = c[0], c1 = c[1];
					for (int n = 0; n < N.size(); n++) {
						if (c0 == N[n])  lc = n;
						if (c1 == N[n])  uc = n;
					}
					for (int i = lc; i <= uc; i++) {
						W1.push_back(W[i]);
					}
				}
			}
			else {
				Split_String(domainparser[t], b, "~");
				string b0, b1;
				b0 = b[0], b1 = b[1];
				for (int n = 0; n < N.size(); n++) {
					if (b0 == N[n])lb = n;
					if (b1 == N[n])ub = n;
				}
				for (int i = lb; i <= ub; i++) {
					W1.push_back(W[i]);
				}
			}
		}
	}
	else {
		vector<string> b, split_frag;
		if (string::size_type id = s.find(",") != string::npos) {
			Split_String(s, b, ",");
			for (unsigned int k = 0; k < b.size(); k++) {
				vector<string> c;
				Split_String(b[k], c, "~");
				string c0, c1;
				c0 = c[0], c1 = c[1];
				for (int n = 0; n < N.size(); n++) {
					if (c0 == N[n])lc = n;
					if (c1 == N[n])uc = n;
				}
				for (int i = lc; i <= uc; i++) {
					W1.push_back(W[i]);
				}
			}
		}
		else {
			Split_String(s, b, "~");
			string b0, b1;
			b0 = b[0], b1 = b[1];
			for (int n = 0; n < N.size(); n++) {
				if (b0 == N[n]) lb = n;
				if (b1 == N[n]) ub = n;
			}
			for (int i = lb; i <= ub; i++) {
				W1.push_back(W[i]);
			}
		}
	}
	return W1;
}
/*
Transform the vetor to string
para:
	vector<string> N -- the initial atom_num
	vector<int> v--the miss element Numof N
	int c1-- the start num of string
	int c2 --the end num of string
return:
	string s -- the string that be successfully spliced
*/
string vector_to_string(vector<string> N, vector<int> v, int c1, int c2)
{
	string s;
	if (v.size() != 0) {
		vector<string> frag;
		string sf = N[c1] + "~" + N[v[0] - 1];
		frag.push_back(sf);
		int k = 0;
		if (v.size() > 1) {
			do {
				int b = k + 1;
				int profix = int(atof(N[v[k]].substr(0, 4).c_str()));
				int suffix = int(atof(N[v[b]].substr(0, 4).c_str()));
				if (suffix != profix + 1) {
					int p = v[k] + 1, s = v[b] - 1;
					string sm = N[p] + "~" + N[b];
					frag.push_back(sm);
				}
				k++;
			} while (k < v.size() - 1);
		}
		int last = v.back() + 1;
		string se = N[last] + "~" + N[c2 - 1];
		frag.push_back(se);
		for (int i = 0; i < frag.size() - 1; i++) {
			s += (frag[i] + ",");
		}
		s += frag.back();
	}
	else {
		vector<int> Dis_continue;
		vector<string> frag;
		for (int t = c1; t < c2 - 1; t++) {
			int profix = int(atof(N[t].substr(0, 4).c_str()));
			int suffix = int(atof(N[t + 1].substr(0, 4).c_str()));
			if (suffix == profix + 1 || suffix == profix) {
				continue;
			}
			else if (suffix != profix + 1) {
				Dis_continue.push_back(t);
			}
		}
		if (Dis_continue.size() != 0) {
			string sf = N[c1] + "~" + N[Dis_continue[0]];
			frag.push_back(sf);
			int i = 0;
			if (Dis_continue.size() > 1) {
				do {
					string sm = N[Dis_continue[i] + 1] + "~" + N[Dis_continue[i + 1]];
					frag.push_back(sm);
					i++;
				} while (i < Dis_continue.size() - 1);
			}
			string se = N[Dis_continue.back() + 1] + "~" + N[c2 - 1];
			frag.push_back(se);
			for (int i = 0; i < frag.size() - 1; i++) {
				s += (frag[i] + ",");
			}
			s += frag.back();
		}
		else {
			s = N[c1] + "~" + N[c2 - 1];
		}
	}
	return s;
}
/*
    generate a 2d array
*/
template <class A> void NewArray(A*** array, int Narray1, int Narray2)
{
	*array = new A * [Narray1];
	for (int i = 0; i < Narray1; i++) *(*array + i) = new A[Narray2];
};

/*
Calculate the distance between two atoms
*/
double distance(double x[3], double y[3])
{
	double d1 = x[0] - y[0];
	double d2 = x[1] - y[1];
	double d3 = x[2] - y[2];
	double d = d1 * d1 + d2 * d2 + d3 * d3;
	return pow(d, 1.0 / 2.0);
}

/*
Merge two vectors
*/
void mergeVector(vector<int> &a, vector<int> &b, vector<int> &c)
{
    merge(a.begin(), a.end(), b.begin(),b.end(), c.begin());
}


/*
Calculate the distance matrix
parameters:
		double **a:2D array
		vector<double> x, vector<double> y, vector<double> z: the vector of X,Y,Z
		int n:length
return:
		double **a
*/
void Cal_Dis_Matrix(double** a, vector<double> x, vector<double> y, vector<double> z, int n)
{
	int i, j = 0;
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			double atomI[3] = { x[i],y[i],z[i] };
			double atomJ[3] = { x[j],y[j],z[j] };
			double dist = distance(atomI, atomJ);
			a[i][j] = dist;
		}
	}
}
/*
Calculate the contact matrix
parameters:
		double **a:2D array
		vector<double> x, vector<double> y, vector<double> z: the vector of X,Y,Z position
		int n:n*n 2D array
return:
		a contact matrix
*/

double Cal_Con_Matrix01(double** a, vector<double> x, vector<double> y, vector<double> z) {
	double d0 = 8.0;
	int length = int(x.size());
	for (int i = 0; i < length; i++) {

		for (int j = 0; j <= i; j++) {

			double atomI[3] = { x[i],y[i],z[i] };
			double atomJ[3] = { x[j],y[j],z[j] };
			double dist = distance(atomI, atomJ);
			if (dist < d0) a[i][j] = 1.0;
			else a[i][j] = 0.0;
			a[j][i] = a[i][j];
		}
	}
	
	return **a;
}
/*
Calculate the contact matrix
parameters:
		double **a:2D array
		vector<double> x, vector<double> y, vector<double> z: the vector of X,Y,Z position
		int n:n*n 2D array
return:
		a contact matrix
*/
double Cal_Con_Matrix(double** a,double** DistanceMatrix,int length) {
	double d0 = 8.0, delta = 1.5;
	for (int i = 0; i < length; i++) {
		for (int j = i; j <length; j++) {
			if (DistanceMatrix[i][j] < 20) {
				double diffDist = (DistanceMatrix[i][j] - d0) / delta;
				double expDiffDist = 1 + exp(diffDist);
				a[i][j] = 1.0 / (expDiffDist);
				
			}
			else {
				a[i][j] = 0;
			}
			a[j][i] = a[i][j];
		}
	}
}


/*
Calculate the contact matrix use contact pro
parameters:
		double **a:2D array
		double **DistanceMatrix: the inter-residue distance
		int length:n*n 2D array
		double** ProMatrix: the inter-residue contact pro
return:
		a contact matrix
*/
double Cal_Con_Matrix_seq(double** a,double** DistanceMatrix,int length,double** ProMatrix) {
	double d0 = 8.0, delta = 1.5;
	for (int i = 0; i < length; i++) {
		for (int j = i; j <length; j++) {
			if (DistanceMatrix[i][j] < 16) {
				double diffDist = (DistanceMatrix[i][j] - d0) / delta;
				double expDiffDist = 1 + exp(diffDist);
				a[i][j] = ProMatrix[i][j] / (expDiffDist);
				
			}
			else {
				a[i][j] = 0;
			}
			a[j][i] = a[i][j];
		}
	}
}

/*
Compute the avedensity of a single domain
i.e. the intral nnc
para:
	ConMatrix
	vector<int> W---the Vector corresponding to the coordinate
return:
	double AveDensity
*/
double Ave_Density(double** ConMatrix, vector<int> W)
{
	int length = int(W.size());
	double density = 0.0, beta = 0.95;
	for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < i-2; j++)
		{
			density += ConMatrix[W[i]][W[j]];
		}
	}
	double AveDensity = density /pow(length,beta);
	return AveDensity;
}
/*
 The inter-nnc of two domains
 para:
	double** a--contact matrix
	vector<int> X1--the coordinate of domain1
	vector<int> X2--the coordinate of domain2
return:
	double nnc
*/
double Inter_nnc(double** a, vector<int> X1, vector<int> X2)
{
	double nc = 0.0, alpha = 0.43;
	for (unsigned int i = 0; i < X1.size(); i++) {
		for (unsigned int j = 0; j < X2.size(); j++) {
			nc += a[X1[i]][X2[j]];
		}
	}
	double nnc = nc / (pow(X1.size(), alpha) * pow(X2.size(), alpha));
	return nnc;
}
/*
Cut the chain into sub-domain
para:
	double **ConMatrix
	vector<int> W---the Vector corresponding to the coordinate
	vector<string> N---the Vector corresponding to the Residue serial number
	vector<int> V-- the second vector
return:
	string domain;
*/
string Cut_Domain_sec(double** ConMatrix, vector<int> W, vector<string> N,vector<int> V)
{
	int length = int(W.size());
	double alpha = 0.43;
	double threshold = Ave_Density(ConMatrix, W)/2;
	//Step1:Cut a Domain into two domians
	//Step1.1:Approach1 One_cut
	double MINnnc1 = 1000;
	int cut1_num = 0;
	string cut1;
	if (length > 60) {
		for (int k = 30; k < length - 30; k++) {
			if (V[k] == 1 ) {
				double nc = 0;
				for (int i = 0; i < k; i++) {
					for (int j = k; j < length; j++) {
						nc += ConMatrix[W[i]][W[j]];
					}
				}
				double nnc;
				nnc = nc / (pow(k, alpha) * pow(length - k, alpha));
				if (nnc < MINnnc1) {
					MINnnc1 = nnc;
					cut1 = N[k];
					cut1_num = k;
				}
			}
		}
	}
	//step1.2:Approach2: two_cut
	string cut21, cut22;
	int cut21_num, cut22_num;
	double MINnnc2 = 10000;
	if (length > 60) {
		for (int k1 = 15; k1 < length - 15; k1++) {
			if (V[k1] == 1 ) {
				for (int k2 = k1 + 35; k2 < length - 15; k2++) {
					if (V[k2] == 1 ) {
						if (ConMatrix[W[k1]][W[k2]] > 0.5) {
							double nc = 0.0;
							for (int i = k1 + 1; i < k2; i++) {
								for (int j = 0; j <= k1; j++) {
									nc += ConMatrix[W[i]][W[j]];
								}
								for (int j = k2; j < length; j++) {
									nc += ConMatrix[W[i]][W[j]];
								}
							}
							double nnc;
							nnc = nc / (pow((k2 - k1 - 1), alpha) * pow((length - k2 + k1 + 1), alpha));
							if (MINnnc2 > nnc) {
								MINnnc2 = nnc;
								cut21 = N[k1];
								cut22 = N[k2];
								cut21_num = k1;
								cut22_num = k2;
							}
						}
					}
				}
			}
		}
	}
	//Step2:Ouput Domain
	string domain;
	if (min(MINnnc1, MINnnc2) < threshold) {
		if (MINnnc1 < MINnnc2) {
			//one_cut
			vector<int> miss_N1_num, miss_N2_num;
			string s1, s2;
			for (int i = 0; i < cut1_num; i++) {
				if (is_string_in_vector(N, N[i]) != 1) {
					miss_N1_num.push_back(i);
				}
			}
			for (int i = cut1_num; i < length; i++) {
				if (is_string_in_vector(N, N[i]) != 1) {
					miss_N2_num.push_back(i);
				}
			}
			s1 = vector_to_string(N, miss_N1_num, 0, cut1_num);
			s2 = vector_to_string(N, miss_N2_num, cut1_num, length);
			domain = s1 + ";" + s2;
		}
		else {
			//two_cut
			vector<int> miss_N1_num, miss_N2_num, miss_N3_num;
			string s1, s2, s3;
			for (int i = cut21_num + 1; i < cut22_num; i++) {
				if (is_string_in_vector(N, N[i]) != 1) {
					miss_N1_num.push_back(i);
				}
			}
			for (int i = 0; i < cut21_num + 1; i++) {
				if (is_string_in_vector(N, N[i]) != 1) {
					miss_N2_num.push_back(i);
				}
			}
			for (int i = cut22_num; i < length; i++) {
				if (is_string_in_vector(N, N[i]) != 1) {
					miss_N3_num.push_back(i);
				}
			}
			s1 = vector_to_string(N, miss_N1_num, cut21_num + 1, cut22_num);
			s2 = vector_to_string(N, miss_N2_num, 0, cut21_num + 1);
			s3 = vector_to_string(N, miss_N3_num, cut22_num, length);
			domain = s1 + ";" + s2 + "," + s3;
		}
	}
	else {
		vector<int> miss_N1;
		string s;
		for (int i = 0; i < length - 1; i++) {
			if (is_string_in_vector(N, N[i]) != 1) {
				miss_N1.push_back(i);
			}
		}
		s = vector_to_string(N, miss_N1, 0, length);
		domain = s;
	}
	//cout << domain << endl;
	return domain;
}
/*

*/
string Cut_Domain(double** ConMatrix, vector<int> W, vector<string> N)
{
	int length = int(W.size());
	double alpha = 0.43;
	double threshold = Ave_Density(ConMatrix, W) / 2;
	double MINnnc1 = 1000;
	int cut1_num = 0;
	if (length > 60) {
		for (int k = 30; k < length - 30; k++) {
			double nc = 0;
			for (int i = 0; i < k; i++) {
				for (int j = k; j < length; j++) {
					nc += ConMatrix[W[i]][W[j]];
				}
			}
			double nnc;
			nnc = nc / (pow(k, alpha) * pow(length - k, alpha));
			if (nnc < MINnnc1) {
				MINnnc1 = nnc;
				cut1_num = k;
			}
		}
	}
	//cout << "MINnnc1:" << MINnnc1 << "    "<<"cut1_num:"<< cut1_num<<endl;
	//step1.2:Approach2: two_cut
	string cut21, cut22;
	int cut21_num, cut22_num;
	double MINnnc2 = 10000;
	if (length > 60) {
		for (int k1 = 15; k1 < length - 15; k1++) {
			for (int k2 = k1 + 35; k2 < length - 15; k2++) {
				if (ConMatrix[W[k1]][W[k2]]>0.5) {
					double nc = 0.0;
					for (int i = k1 + 1; i < k2; i++) {
						for (int j = 0; j <= k1; j++) {
							nc += ConMatrix[W[i]][W[j]];
						}
						for (int j = k2; j < length; j++) {
							nc += ConMatrix[W[i]][W[j]];
						}
					}
					double nnc;
					nnc = nc / (pow((k2 - k1 - 1), alpha) * pow((length - k2 + k1 + 1), alpha));
					if (MINnnc2 > nnc) {
						MINnnc2 = nnc;
						cut21_num = k1;
						cut22_num = k2;
					}
				}
			}
		}
	}
	//cout << "MINnnc2:" << MINnnc2 << "    "<<cut21_num <<","<< cut22_num <<endl;
	//cout << threshold << endl;
	//Step2:Ouput Domain
	string domain;
	if (min(MINnnc1, MINnnc2) < threshold) {
		if (MINnnc1 < MINnnc2) {
			//one_cut
			vector<int> miss_N1_num, miss_N2_num;
			string s1, s2;
			for (int i = 0; i < cut1_num; i++) {
				if (is_string_in_vector(N, N[i]) != 1) {
					miss_N1_num.push_back(i);
				}
			}
			for (int i = cut1_num; i < length; i++) {
				if (is_string_in_vector(N, N[i]) != 1) {
					miss_N2_num.push_back(i);
				}
			}
			s1 = vector_to_string(N, miss_N1_num, 0, cut1_num);
			s2 = vector_to_string(N, miss_N2_num, cut1_num, length);
			domain = s1 + ";" + s2;
		}
		else {
			//two_cut
			vector<int> miss_N1_num, miss_N2_num, miss_N3_num;
			string s1, s2, s3;
			for (int i = cut21_num + 1; i < cut22_num; i++) {
				if (is_string_in_vector(N, N[i]) != 1) {
					miss_N1_num.push_back(i);
				}
			}
			for (int i = 0; i < cut21_num + 1; i++) {
				if (is_string_in_vector(N, N[i]) != 1) {
					miss_N2_num.push_back(i);
				}
			}
			for (int i = cut22_num; i < length; i++) {
				if (is_string_in_vector(N, N[i]) != 1) {
					miss_N3_num.push_back(i);
				}
			}
			s1 = vector_to_string(N, miss_N1_num, cut21_num + 1, cut22_num);
			s2 = vector_to_string(N, miss_N2_num, 0, cut21_num + 1);
			s3 = vector_to_string(N, miss_N3_num, cut22_num, length);
			domain = s1 + ";" + s2 + "," + s3;
		}
	}
	else {
		vector<int> miss_N1;
		string s;
		for (int i = 0; i < length - 1; i++) {
			if (is_string_in_vector(N, N[i]) != 1) {
				miss_N1.push_back(i);
			}
		}
		s = vector_to_string(N, miss_N1, 0, length);
		domain = s;
	}
	//cout << domain << endl;
	return domain;
}
vector<string> Combine_Domain(vector<string> initialFrag,vector<int> W,vector<string> N,double **ContactMatrix) {
	vector<string> domainsOut;
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
			vector<int> W1;
			W1 = string_to_vector_int(initialFrag[i], W, N);
			//THE Intra nnc i.e the domain's average density
			double W1nnc = Ave_Density(ContactMatrix, W1) /2;
			for (unsigned int j = 0; j < i; j++) {
				//The second sub_domain
				vector<int> W2;
				cout << initialFrag[i] << "/"<< initialFrag[j]<<endl;
				W2 = string_to_vector_int(initialFrag[j], W, N);
				//THE Intra nnc i.e the domain's average density
				double W2nnc = Ave_Density(ContactMatrix, W2) /2;
				//THE Inter nnc
				double nnc_W1_W2 = Inter_nnc(ContactMatrix, W1, W2);
				double MINnnc = min(W1nnc, W2nnc);
				//Judge the two domain whether can be merged
				if (nnc_W1_W2 >= MINnnc) {
					// Judge the split condition 
					// step1: concat the fragment
					vector<int> W1W2;
					// W1W2 = string_to_vector_int(initialFrag[i] + "," + initialFrag[j],W,N);
					W1W2.resize(W1.size()+W2.size());
					mergeVector(W1,W2,W1W2);
					// step2: the W1W2nnc
					double W1W2nnc = Ave_Density(ContactMatrix, W1W2) /2;
					if(nnc_W1_W2 >= W1W2nnc) //no statisfied the cut condition
					{
						double dif_nnc = nnc_W1_W2 - MINnnc; //Amount of change
						if (MaxNnc < dif_nnc) 
						{
							MaxNnc = dif_nnc;
							combineFrag = initialFrag[i] + ";" + initialFrag[j];
						}
					}
				}
			}
		}
		//cout << "CombineFrag:" << combineFrag << endl;
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
		else{
			domainsOut.swap(initialFrag);
			q = 0;
		}
	} while (q > 0);
	return domainsOut;
}
