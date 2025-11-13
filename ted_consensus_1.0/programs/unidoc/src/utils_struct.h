#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>
#include <string.h>
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
Judege the element whether in vector
para:
	vector<int> v
	int e
*/
bool is_element_in_vector(vector<int> v, int e) {
	vector<int>::iterator it;
	it = find(v.begin(), v.end(), e);
	if (it != v.end()) {
		return true;
	}
	else {
		return false;
	}
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
	vector<int> W--the vector of initial domain
	vector<int> N--the Residue serial number of  initial domain
return:
	vector<int> W1 ---the vector of string
*/
vector<int> string_to_vector(const string s, const vector<int> W, const vector<int> N)
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
					int c0 = int(atof(c[0].c_str())), c1 = int(atof(c[1].c_str()));
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
				int b0, b1;
				b0 = int(atof(b[0].c_str())), b1 = int(atof(b[1].c_str()));
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
				int c0, c1;
				c0 = int(atof(c[0].c_str())), c1 = int(atof(c[1].c_str()));
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
			int b0, b1;
			b0 = int(atof(b[0].c_str())), b1 =int(atof(b[1].c_str()));
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
	vector<int> v--the miss element of N
	int c1-- the start of string
	int c2 --the end of string
return:
	string s -- the string that be successfully spliced
*/
string vector_to_string(vector<int> v, int c1, int c2)
{
	vector<string> frag;
	string s;
	if (v.size() != 0) {
		string sf = to_string(c1) + "~" + to_string(v[0] - 1);
		frag.push_back(sf);
		int k = 0;
		if (v.size() > 1) {
			do {
				if (v[k + 1] != v[k] + 1) {
					string sm = to_string(v[k] + 1) + "~" + to_string(v[k + 1] - 1);
					frag.push_back(sm);
				}
				k++;
			} while (k < v.size() - 1);
		}
		string se = to_string(v.back() + 1) + "~" + to_string(c2 - 1);
		frag.push_back(se);
		for (int i = 0; i < frag.size() - 1; i++) {
			s += (frag[i] + ",");
		}
		s += frag.back();
	}
	else {
		s = to_string(c1) + "~" + to_string(c2 - 1);
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
Calculate the distance matrix
parameters:
		double **a:2D array
		vector<double> x, vector<double> y, vector<double> z: the vector of X,Y,Z
		int n:n*n 2D array
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

double Cal_Con_Matrix(double** a, vector<double> x, vector<double> y, vector<double> z) {
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
	double density = 0.0, alpha = 0.43;
	for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < i - 2; j++)
		{
			density += ConMatrix[W[i]][W[j]];
		}
	}
	double AveDensity = density / (length);
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
	vector<int> N---the Vector corresponding to the Residue serial number
return:
	string domain;
*/
string Cut_Domain_sec(double** ConMatrix, vector<int> W, vector<int> N,vector<int> V)
{
	int length = int(W.size());
	double alpha = 0.43;
	double threshold = Ave_Density(ConMatrix, W) / 2;
	//Step1:Cut a Domain into two domians
	//Step1.1:Approach1 One_cut
	double MINnnc1 = 1000;
	int cut1 = 0;
	if (length > 60) {
		for (int k = 30; k < length - 30; k++) {
			double nc = 0;
			if (V[k] == 1) {
				//The second structure is Coil and Turn will be cut
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
				}
			}
		}
	}
	//step1.2:Approach2: two_cut
	int cut21, cut22;
	double MINnnc2 = 10000;
	if (length > 60) {
		for (int k1 = 15; k1 < length - 15; k1++) {
			if (V[k1] == 1  ) {
				for (int k2 = k1 + 35; k2 < length - 15; k2++) {
					if (V[k2] == 1 ) {
						if (ConMatrix[W[k1]][W[k2]] > 0.5) {
							double nc = 0.0;
							for (int i = k1 + 1; i < k2; i++) {
								for (int j = 0; j < k1; j++) {
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
							}
						}
					}
				}
			}
		}
	}
	//Step2:Ouput Domain
	string domain;
	if (min(MINnnc1, MINnnc2) <= threshold) {
		if (MINnnc1 < MINnnc2) {
			//one_cut
			vector<int> miss_N1, miss_N2;
			string s1, s2;
			for (int i = N[0]; i < cut1; i++) {
				if (is_element_in_vector(N, i) != 1) {
					miss_N1.push_back(i);
				}
			}
			for (int i = cut1; i <= N[length - 1]; i++) {
				if (is_element_in_vector(N, i) != 1) {
					miss_N2.push_back(i);
				}
			}
			s1 = vector_to_string(miss_N1, N[0], cut1);
			s2 = vector_to_string(miss_N2, cut1, N.back() + 1);
			domain = s1 + ";" + s2;
		}
		else {
			//two_cut
			vector<int> miss_N1, miss_N2, miss_N3;
			string s1, s2, s3;
			for (int i = cut21 + 1; i < cut22; i++) {
				if (is_element_in_vector(N, i) != 1) {
					miss_N1.push_back(i);
				}
			}
			for (int i = N[0]; i < cut21 + 1; i++) {
				if (is_element_in_vector(N, i) != 1) {
					miss_N2.push_back(i);
				}
			}
			for (int i = cut22; i <= N[length - 1]; i++) {
				if (is_element_in_vector(N, i) != 1) {
					miss_N3.push_back(i);
				}
			}
			s1 = vector_to_string(miss_N1, cut21 + 1, cut22);
			s2 = vector_to_string(miss_N2, N[0], cut21 + 1);
			s3 = vector_to_string(miss_N3, cut22, N.back() + 1);
			domain = s1 + ";" + s2 + "," + s3;
		}
	}
	else {
		vector<int> miss_N1;
		string s;
		for (int i = N[0]; i < N[length - 1]; i++) {
			if (is_element_in_vector(N, i) != 1) {
				miss_N1.push_back(i);
			}
		}
		s = vector_to_string(miss_N1, N[0], N.back() + 1);
		domain = s;
	}
	return domain;
}
string Cut_Domain(double** ConMatrix, vector<int> W, vector<int> N)
{
	int length = int(W.size());
	double alpha = 0.43;
	double threshold = Ave_Density(ConMatrix, W) / 2;
	//Step1:Cut a Domain into two domians
	//Step1.1:Approach1 One_cut
	double MINnnc1 = 1000;
	int cut1 = 0;
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
				cut1 = N[k];
			}
		}
	}
	//step1.2:Approach2: two_cut
	int cut21, cut22;
	double MINnnc2 = 10000;
	if (length > 60) {
		for (int k1 = 15; k1 < length - 15; k1++) {
			for (int k2 = k1 + 35; k2 < length - 15; k2++) {
				if (ConMatrix[W[k1]][W[k2]] == 1) {
					double nc = 0.0;
					for (int i = k1 + 1; i < k2; i++) {
						for (int j = 0; j < k1; j++) {
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
			vector<int> miss_N1, miss_N2;
			string s1, s2;
			for (int i = N[0]; i < cut1; i++) {
				if (is_element_in_vector(N, i) != 1) {
					miss_N1.push_back(i);
				}
			}
			for (int i = cut1; i <= N[length - 1]; i++) {
				if (is_element_in_vector(N, i) != 1) {
					miss_N2.push_back(i);
				}
			}
			s1 = vector_to_string(miss_N1, N[0], cut1);
			s2 = vector_to_string(miss_N2, cut1, N.back() + 1);
			domain = s1 + ";" + s2;
		}
		else {
			//two_cut
			vector<int> miss_N1, miss_N2, miss_N3;
			string s1, s2, s3;
			for (int i = cut21 + 1; i < cut22; i++) {
				if (is_element_in_vector(N, i) != 1) {
					miss_N1.push_back(i);
				}
			}
			for (int i = N[0]; i < cut21 + 1; i++) {
				if (is_element_in_vector(N, i) != 1) {
					miss_N2.push_back(i);
				}
			}
			for (int i = cut22; i <= N[length - 1]; i++) {
				if (is_element_in_vector(N, i) != 1) {
					miss_N3.push_back(i);
				}
			}
			s1 = vector_to_string(miss_N1, cut21 + 1, cut22);
			s2 = vector_to_string(miss_N2, N[0], cut21 + 1);
			s3 = vector_to_string(miss_N3, cut22, N.back() + 1);
			domain = s1 + ";" + s2 + "," + s3;
		}
	}
	else {
		vector<int> miss_N1;
		string s;
		for (int i = N[0]; i < N[length - 1]; i++) {
			if (is_element_in_vector(N, i) != 1) {
				miss_N1.push_back(i);
			}
		}
		s = vector_to_string(miss_N1, N[0], N.back() + 1);
		domain = s;
	}
	return domain;
}
