// 02_one_loop_renormalization.cc
//
// Color factors of the one-loop QCD diagrams in the journal notebook
//   msc-qcd-journal/notes-QFT/RenormalizationQCD.nb
// That notebook does the Lorentz/loop part with FeynCalc. Here we redo only
// the color part with ColorFull. Run with:
//   make run EX=02_one_loop_renormalization
//
// Diagrams covered (same as diagGluonSE, diagGhostSE, diagQuarkSE,
// diagQuarkGluonVertex in the notebook):
//   A) gluon self-energy, quark bubble          Tr(t^a t^b)      = TR delta^ab
//   B) gluon self-energy, gluon/ghost bubble     f^acd f^bcd      = CA delta^ab
//   C) quark self-energy                         t^c t^c          = CF 1
//   D) quark-gluon vertex, abelian diagram       t^c t^a t^c      = (CF - CA/2) t^a
//   E) quark-gluon vertex, non-abelian diagram   f^abc t^b t^c    = (i/2) CA t^a
// and the one-loop beta-function coefficient built from them,
//   beta0 = 11/3 CA - 4/3 TR nf.
//
// How f^abc is written in ColorFull (trace basis only knows t's):
//   Tr(t^a t^b t^c) - Tr(t^a t^c t^b) = Tr(t^a [t^b,t^c]) = i TR f^abc
// so with F(a,b,c) := (a,b,c) - (a,c,b) we have f^abc = -(i/TR) F(a,b,c).
//
// A repeated label is summed over (contracted), like a repeated index.
//
// How a number is read off: if X = c * B for a known structure B, then
//   <B|X> = c <B|B>, so c = <B|X> / <B|B>.
// Normalization: Tr(t^a t^b) = TR delta^ab and [t^a,t^b] = i f^abc t^c give
// CA = 2 TR Nc, which is the familiar CA = Nc only for TR = 1/2 (the default).
//
// We do this at Nc = 3 and Nc = 5. CA and CF = TR(Nc^2-1)/Nc change
// differently with Nc, so two values of Nc tell them apart.

#include <iostream>
#include <iomanip>
#include "ColorFull/Core/Col_functions.h"
#include "ColorFull/Core/Col_amp.h"

using namespace ColorFull;
using namespace std;

// F(a,b,c) = Tr(t^a t^b t^c) - Tr(t^a t^c t^b) = i TR f^abc
Col_amp F(int a, int b, int c)
{
	auto s = [](int x, int y, int z)
	{
		return "[(" + to_string(x) + "," + to_string(y) + "," + to_string(z) + ")]";
	};
	return Col_amp(s(a, b, c)) - Col_amp(s(a, c, b));
}

// c such that X = c * B (numerically, at the Nc currently set in cf)
double coefficient(const Col_functions &cf, const Col_amp &B, const Col_amp &X)
{
	return cf.double_num(cf.scalar_product(B, X)) / cf.double_num(cf.scalar_product(B, B));
}

int main()
{
	Col_functions cf; // TR = 1/2 by default

	// Two-gluon structure (a,b) = Tr(t^a t^b) = TR delta^ab.  Gluons: a=1, b=2.
	Col_amp gluon_ab("[(1,2)]");
	// Quark structures.  Quark = 1, antiquark = 2, external gluon a = 3.
	Col_amp delta_q("[{1,2}]"); // delta_ij
	Col_amp ta_q("[{1,3,2}]");  // (t^a)_ij

	// ---- A) quark bubble: t^a_ij t^b_ji = Tr(t^a t^b) -----------------------
	Col_amp quark_bubble("[(1,2)]");

	// ---- B) gluon bubble: f^acd f^bcd, with c=3, d=4 summed ----------------
	// f^acd f^bcd = (-i/TR)^2 F(1,3,4) F(2,3,4) = -(1/TR^2) F F
	Col_amp FF = F(1, 3, 4) * F(2, 3, 4);

	// ---- C) quark self-energy: (t^c t^c)_ij, c = 4 summed -----------------
	Col_amp quark_se("[{1,4,4,2}]");

	// ---- D) abelian vertex: (t^c t^a t^c)_ij, c = 4 summed ----------------
	Col_amp vertex_abelian("[{1,4,3,4,2}]");

	// ---- E) non-abelian vertex: f^abc (t^b t^c)_ij, b=4, c=5 summed ------
	// f^abc t^b t^c = -(i/TR) F(3,4,5) {1,4,5,2}
	Col_amp vertex_nonabelian = F(3, 4, 5) * Col_amp("[{1,4,5,2}]");

	// Exact polynomials. (ColorFull's contract_all_gluons needs every label
	// to appear twice, so we close the open a,b by contracting with (a,b).)
	// If f^acd f^bcd = CA delta^ab then
	//   <(a,b)| F F> = -TR^2 * CA * TR (Nc^2-1) = -2 TR^4 Nc (Nc^2-1).
	// ColorFull writes everything in TR, Nc, CF; use Nc CF = TR(Nc^2-1) to compare.
	cout << "Exact ColorFull polynomials (all labels contracted):\n"
	     << "  <(a,b)|(a,b)>        = " << cf.scalar_product(gluon_ab, gluon_ab)
	     << "   [TR^2 (Nc^2-1)]\n"
	     << "  <(a,b)|F(acd)F(bcd)> = " << cf.scalar_product(gluon_ab, FF)
	     << "   [-TR^2 CA TR(Nc^2-1) = -2 TR^4 Nc (Nc^2-1)]\n"
	     << "  <{q,qb}|t^c t^c>     = " << cf.scalar_product(delta_q, quark_se)
	     << "   [Nc CF]\n"
	     << "  <t^a|t^c t^a t^c>    = " << cf.scalar_product(ta_q, vertex_abelian)
	     << "   [TR(Nc^2-1)(CF-CA/2) = -TR^2 (Nc^2-1)/Nc]\n\n";

	cout << fixed << setprecision(4);
	for (double Nc : {3.0, 5.0})
	{
		cf.set_Nc(Nc);
		double TR = cf.get_TR();
		double CA = 2 * TR * Nc; // = Nc for TR = 1/2
		double CF = TR * (Nc * Nc - 1) / Nc;
		// X = c delta^ab = (c/TR) (1,2), so c = TR * coefficient(...)
		double A = TR * coefficient(cf, gluon_ab, quark_bubble);
		double B = -1.0 / (TR * TR) * TR * coefficient(cf, gluon_ab, FF);
		double C = coefficient(cf, delta_q, quark_se);
		double D = coefficient(cf, ta_q, vertex_abelian);
		// E is i * E_im; the -(i/TR) prefactor is applied here
		double E_im = -1.0 / TR * coefficient(cf, ta_q, vertex_nonabelian);

		cout << "Nc = " << Nc << "  (TR = " << TR << ", CA = " << CA
		     << ", CF = " << CF << ")\n";
		cout << "  A) quark bubble     Tr(t^a t^b)  = " << setw(8) << A
		     << " delta^ab   expect TR        = " << TR << "\n";
		cout << "  B) gluon bubble     f^acd f^bcd  = " << setw(8) << B
		     << " delta^ab   expect CA        = " << CA << "\n";
		cout << "  C) quark self-en.   t^c t^c      = " << setw(8) << C
		     << " 1          expect CF        = " << CF << "\n";
		cout << "  D) abelian vertex   t^c t^a t^c  = " << setw(8) << D
		     << " t^a        expect CF-CA/2   = " << CF - CA / 2 << "\n";
		cout << "  E) non-ab. vertex   f^abc t^b t^c= " << setw(7) << E_im
		     << "i t^a        expect (i/2)CA   = " << CA / 2 << "i\n";
		for (int nf : {3, 5, 6})
		{
			double beta0 = 11.0 / 3.0 * B - 4.0 / 3.0 * A * nf;
			cout << "  beta0(nf=" << nf << ") = 11/3 CA - 4/3 TR nf = " << beta0 << "\n";
		}
		cout << "\n";
	}
	return 0;
}
