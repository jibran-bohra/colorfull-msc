// 01_basic_color_factors.cc
//
// First ColorFull program: textbook SU(Nc) color factors that you can
// verify by hand. Run with:   make run EX=01_basic_color_factors
//
// ColorFull notation (trace basis):
//   {q, g1, g2, ..., qbar}  an open quark line: (t^g1 t^g2 ...)_{q qbar}
//   (g1, g2, ..., gn)       a closed gluon loop: Tr(t^g1 t^g2 ... t^gn)
//   [ ... ]                 a product of such objects (one "color structure")
//
// Results are polynomials in Nc, TR (Tr(t^a t^b) = TR delta^ab) and CF.
// Numerical values use Nc = 3, TR = 1/2, CF = 4/3 by default.

#include <iostream>
#include "ColorFull/Core/Col_functions.h"
#include "ColorFull/Core/Col_amp.h"
#include "ColorFull/Core/Trace_basis.h"

using namespace ColorFull;
using namespace std;

int main() {
	Col_functions cf;  // the toolbox: scalar products, gluon exchange, ...

	// ---------------------------------------------------------------
	// 1) q qbar: the color singlet delta_{q qbar}.
	//    |delta|^2 = sum_{ij} delta_ij delta_ij = Nc
	Col_amp qq("[{1,2}]");
	cout << "1) <qq|qq>             = " << cf.scalar_product(qq, qq)
	     << "    (expect Nc)" << endl;

	// ---------------------------------------------------------------
	// 2) q qbar g: the amplitude (t^a)_{ij}.
	//    sum |t^a_ij|^2 = Tr(t^a t^a) = TR (Nc^2 - 1) = Nc CF
	Col_amp qqg("[{1,3,2}]");
	Polynomial p2 = cf.scalar_product(qqg, qqg);
	cout << "2) <qqg|qqg>           = " << p2
	     << "    (expect TR(Nc^2-1) = 4 at Nc=3)  numeric: "
	     << cf.double_num(p2) << endl;

	// ---------------------------------------------------------------
	// 3) g g: Tr(t^a t^b) = TR delta^ab, so |.|^2 = TR^2 (Nc^2-1)
	Col_amp gg("[(1,2)]");
	Polynomial p3 = cf.scalar_product(gg, gg);
	cout << "3) <gg|gg>             = " << p3
	     << "    numeric: " << cf.double_num(p3) << "  (expect 2)" << endl;

	// ---------------------------------------------------------------
	// 4) Color correlator <M| T_1 . T_2 |M> for a q qbar singlet.
	//    Color conservation T_1 + T_2 = 0 gives T_1.T_2 = -CF, times <M|M> = Nc.
	Polynomial p4 = cf.color_correlator(qq, 1, 2);
	cout << "4) <qq|T1.T2|qq>       = " << p4
	     << "    numeric: " << cf.double_num(p4)
	     << "  (expect -CF*Nc = -4)" << endl;

	// ---------------------------------------------------------------
	// 5) Three gluons: f^{abc} f^{abc} = Nc (Nc^2 - 1) = 24 at Nc = 3.
	//    With i f^{abc} TR = Tr(t^a t^b t^c) - Tr(t^a t^c t^b) we need
	//    |Tr(123) - Tr(132)|^2 = <123|123> + <132|132> - 2 Re<123|132>.
	Col_amp t123("[(1,2,3)]");
	Col_amp t132("[(1,3,2)]");
	Polynomial a = cf.scalar_product(t123, t123);
	Polynomial b = cf.scalar_product(t132, t132);
	Polynomial c = cf.scalar_product(t123, t132);
	double ff = (cf.double_num(a) + cf.double_num(b) - 2 * cf.double_num(c))
	            / (0.5 * 0.5);  // divide by TR^2
	cout << "5) f^abc f^abc         = " << ff << "  (expect 24)" << endl;
	cout << "   <123|123> = " << a << ",  <123|132> = " << c << endl;

	// ---------------------------------------------------------------
	// 6) A whole basis: q qbar g g. ColorFull builds the trace basis
	//    and its scalar-product (Gram) matrix as polynomials in Nc.
	Trace_basis basis(1, 2);  // 1 q qbar pair, 2 gluons
	cout << "\n6) Trace basis for q qbar g g:\n" << basis;
	basis.scalar_product_matrix();
	cout << "Scalar product matrix (exact):\n" << basis.P_spm;
	cout << "Scalar product matrix (Nc=3):\n" << basis.d_spm;

	return 0;
}
