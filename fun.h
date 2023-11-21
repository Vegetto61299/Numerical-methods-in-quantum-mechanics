#pragma once
#include <math.h>
#include <iostream>
#include <armadillo>
#include <fstream>
#include <bitset>
#include <string>
#include <vector>
// pert=h*cos(q*i)*sz
std::complex<double> ii(0, 1);
using namespace std;
using namespace arma;
#define Lmax 16

class hamiltonian {
private:
	int spm(int x, int pos) {
		std::bitset<Lmax> k = x;
		if (k[(pos + 1)%L] == 0)
			return -1;
		else if (k[pos] == 1)
			return -1;
		k[(pos + 1) % L] = 0;
		k[pos] = 1;
		return k.to_ulong();
	}
	int smp(int x, int pos) {
		std::bitset<Lmax> k = x;
		if (k[(pos + 1) % L] == 1)
			return -1;
		else if (k[pos] == 0)
			return -1;
		k[(pos + 1) % L] = 1;
		k[pos] = 0;
		return k.to_ulong();
	}
	double sz(int x, int pos) {
		std::bitset<Lmax> k = x;
		if (k[pos % L] == 1)
			return 0.5;
		else return -0.5;
	}
	void HMake(double q) {
		int s = 0, k = 0;
		for (int i = 0; i < d; i++) {
			for (int j = 0; j < L-1; j++) {//periodyczne warunki brzegowe dla L, otwarte dla L-1
				H(i, i) += delta * sz(i, j) * sz(i, j + 1)+ h * cos(q * j) * macszi(j)(i);
				s = smp(i, j);
				k = spm(i, j);
				if (smp(i, j) > 0)
					H(i, s) += J * 0.5;
				if (spm(i, j) > 0)
					H(i, k) += J * 0.5;
			}
		}
	}




public:
	hamiltonian(int size, double delt, double constj,double h0,double q,int lanczos) {
		delta = delt;
		J = constj;
		L = size;
		h = h0;
		d = pow(2, L);
		H = dmat(d, d, fill::zeros);
		HMake(q);
		if (lanczos==0)
		eig_sym(eigval, eigvec, H);
	}
	double egap() {
		return eigval(1) - eigval(0);
	}

	double pojciep(double T) {
		double p;
		double psum = 0;
		double E2av = 0;
		double Eav2 = 0;
		for (int i = 0; i < d; i++) {
			p = exp(-eigval(i) / T);
			psum += p;
			Eav2 += p * eigval(i);
			E2av += p * eigval(i) * eigval(i);
		}
		E2av /= psum;
		Eav2 /= psum;
		Eav2 *= Eav2;
		return (E2av - Eav2) / T / T / L;
	}

	cx_vec eczas(double t, vec vecp) {
		cx_vec psit(d, fill::zeros);
		for (int i = 0; i < d; i++) {
			psit += dot(eigvec.col(i),vecp)*exp(-ii * t * eigval(i))*eigvec.col(i);
		}
		return psit;
	}

	vec macszi(int i) {//tworzy wektor diagonali Szi
		std::bitset<Lmax> k;
		vec out(d, fill::zeros);
		out -= 0.5;
		for (int j = 0; j < d; j++) {
			k = j;
			if (k[i] == 1) {// da³bym 1 ale 0 bo uk³ad siê odbija tzn zamiast 0011 mamy 1100
				out(j) += 1;
			}
		}
		return out;
	}

	vec magn(cx_vec psi) {// funkcja zwraca œredni¹ magnetyzacjê w stanie psi
		vec praw(d, fill::zeros);//macierz prawdopodobieñstw
		for (int i = 0; i < d; i++) {
			praw(i) = real(psi(i)*conj(psi(i)));
		}
		std::bitset<Lmax> k;
		vec out(L, fill::zeros);
		for (int i = 0; i < d; i++) {
			k = i;
			for (int j = 0; j < L; j++) {
				if (k[j] == 0){// da³bym 1 ale 0 bo uk³ad siê odbija tzn zamiast 0011 mamy 1100
					out(j) += praw(i);
				}
			}
		}
		return out-0.5;
	}

	double delta;
	double J;
	int L;
	int d;
	double h;
	vec eigval;
	mat eigvec;
	dmat H;
	int myNum;
	string myString;
};

void cieplo() {
	hamiltonian* hamk;
	std::fstream plik;
	int N = 100;
	dmat res(6, N, fill::zeros);
	plik.open("data.txt", std::ios::in | std::ios::out);
	for (int j = 0; j < N; j++)
		res(0, j)= j * 4.0 / (N - 1);
	for (int i = 0; i < 7; i++) {
		hamk = new hamiltonian(2 * i + 2, 1, 1, 0, 0,0);
		for (int j = 0; j < N; j++) {
			res(i+1, j) = hamk->pojciep(j * 4.0/(N-1));
		}
	}
	plik << trans(res);
	cout << trans(res);
}
void magnetyzacja() {
	hamiltonian ham(14, 1, 1, 0,0,0);
	vec vecp(ham.d, fill::zeros);
	vecp(pow(2, (ham.L / 2)) - 1) = 1;
	for (int i = 0; i < 200; i++) {
		cout << i * 0.05 << trans(ham.magn(ham.eczas(i*0.05,vecp)));
	}
}



void histogram(double T,int N) {//ten z 89 strony, kod trochê ma³o wydajny
	double beta = 0;
	hamiltonian ham(8, 1, 1, 0, 0, 0);
	double wmax = ham.eigval(ham.d - 1) - ham.eigval(0);
	double bin = wmax / double(N)*2;
	double q,Z=0;
	vec szi[Lmax];
	int w;
	/*vec omegi(N + 1, fill::zeros);
	for (int i = 0; i<N + 1; i++)
		omegi(i) = -wmax+bin*i;
	cout << omegi;*/
	dmat s = dmat(ham.L+1, N+1, fill::zeros);
	std::complex<double> a(0, 0);
	cx_vec szq(ham.d, fill::zeros);
	for (int i = 0; i < ham.d; i++) { //robimy macierz Szq i liczymy Z
		Z += exp(-beta * ham.eigval(i));
	}
	for (int i = 0; i < ham.L; i++)
		szi[i] = ham.macszi(i);
	for (int k = 0; k <= ham.L; k++) {//niepewne ale chyba <=
		q = 4 * acos(0) / double(ham.L)*k;
		szq =szq*0;
		for (int i = 0; i < ham.L; i++) {//robimy macierz Szq
			szq += exp(i*q*ii)* szi[i];
		}
		for (int n = 0; n < ham.d; n++) {
			for (int m = 0; m < ham.d; m++) {//m>n? czy wartoœæ bezwzglêdna
				if (T != 0) {
					w = (ham.eigval(m) - ham.eigval(n) + wmax) / bin;
					a = dot(conj(ham.eigvec.col(n)), szq % ham.eigvec.col(m));
					s(k, w) += exp(-beta * ham.eigval(n)) / Z / bin * real(a * conj(a));
				}
				else if (n == 0) {
					w = (ham.eigval(m) - ham.eigval(n) + wmax) / bin;
					a = dot(conj(ham.eigvec.col(n)), szq % ham.eigvec.col(m));
					s(k, w) += real(a * conj(a)) / bin;
				}
			}
		}
	}
	cout << trans(s);
	}

void perturbacja(int size, double delt, double constj, double h0) {
	hamiltonian h(size, delt, constj, 0, 0,0);
	hamiltonian *hi;
	double tmax = 60,steps=3001;
	double dt = tmax / steps;
	double q,omeg;
	dmat magn(size, steps,fill::zeros);
	cx_dmat res(2, 301, fill::zeros);

	for (int k = 0; k < 2; k++) {
		q = (acos(0) + k * acos(0));
		hi = new hamiltonian(size, delt, constj, h0, q,0);
		for (int t = 0; t * dt < tmax; t++) 
			magn.col(t) = h.magn(h.eczas(t*dt,hi->eigvec.col(0)));
		for (int w = 0; w < 301; w++) {
			omeg = acos(0) * 2.0 / 300.0*w;
			for (int l = 0; l < size; l++) {
				for (int t = 0; t * dt < tmax; t++) {
					res(k, w) += exp(q*l*ii)*dt * exp(omeg * t * dt * ii) * magn(l, t);
				}
			}
			res(k, w) = res(k, w)*conj(res(k, w));
		}
	}
	cout << trans(real(res));
}

void lanczos(int N) {//105 i 106 strona zrób te wykresy

	std::fstream plik;
	plik.open("data.txt", std::ios::in | std::ios::out);

	int L = 11, d = pow(2,L);
	hamiltonian h(L, 1, 1, 0, 0, 0);
	dmat ziemniaczek(N, N, fill::zeros);
	dmat eigvec(N,N, fill::zeros);
	vec fiim1, fii, hfi, tmp2, eigval, result(39,fill::zeros);//fi_i-1
	double ai=0, bi=0;

	fii.randu(d);//krok 0
	fii = (fii * 2 - 1);
	fii = fii / sqrt(dot(fii , fii));

	hfi = h.H * fii;
	ai = dot(fii, hfi);//a0
	tmp2 = hfi - ai * fii;
	ziemniaczek(0, 0) = ai;

	for (int i = 1; i < N; i++) {
		fiim1 = fii;
		ziemniaczek.resize(i + 1,i+1);
		eigvec.clear();
		eigvec.resize(i + 1, i + 1);
		eigval *= 0;
		bi = sqrt(dot(tmp2, tmp2));
		fii = tmp2 / bi;
		hfi = h.H * fii;
		ai = dot(fii, hfi);
		tmp2 = hfi - ai * fii-bi*fiim1;

		ziemniaczek(i, i) = ai;
		ziemniaczek(i-1, i) = bi;
		ziemniaczek(i, i-1) = bi;
		eig_sym(eigval, eigvec, ziemniaczek);
		//plik << trans(eigval);
		result(i-1) = eigval(0);
	}
	cout << (result-eigval(0))/eigval(0);
	//eig_sym(eigval, eigvec, h.H);
	//plik << trans(eigval);
	//cout << trans(eigval);
	//cout << eigval;

}
