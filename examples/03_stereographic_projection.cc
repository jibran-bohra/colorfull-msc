// 03_stereographic_projection.cc
//
// The stereographic projection from the journal notebook
//   msc-qcd-journal/notes-CFT/stereographic.nb
// redone in C++, then joined to ColorFull colour factors. Run with:
//   make run EX=03_stereographic_projection
//
// Part 1, geometry (the notebook's rr[x] and rrd[x]):
//   A point x = (x1,x2) on the transverse plane goes to the unit vector
//     n(x) = ( 2 x1, 2 x2, 1 - x^2 ) / (1 + x^2)
//   i.e. cos(theta) = (1-x^2)/(1+x^2), |x| = tan(theta/2). x = 0 is the north
//   pole, |x| -> infinity is the south pole (the projection point).
//   With a sphere of radius r and the plane at height d (notebook's rrd),
//     w(x) = ( 2r(d+r) x1, 2r(d+r) x2, r((d+r)^2 - x^2) ) / ((d+r)^2 + x^2),
//   which lies on the line from s = (0,0,-r) through (x1,x2,d).
//   We check: |n| = 1, the inverse map, |x| = tan(theta/2), the rrd form,
//   conformality (J^T J = Omega^2 * 1 with Omega = 2/(1+x^2)), and the
//   notebook's key identity
//     1 - n(x).n(y) = 2 (x-y)^2 / ((1+x^2)(1+y^2)).
//
// Part 2, kernels: that identity turns the BMS (jet) kernel on the sphere,
//     W_ab(k) dOmega_k = n_ab / (n_ak n_kb) dOmega_k,   n_ab = 1 - n_a.n_b,
//   into the JIMWLK (CGC) dipole kernel on the plane,
//     K_ab(z) d^2z     = (x_a-x_b)^2 / ((x_a-z)^2 (z-x_b)^2) d^2z.
//   With dOmega = Omega^2 d^2x the (1+x^2) factors all cancel and
//     W_ab dOmega = 2 K_ab d^2z,
//   point by point. We check this ratio at random points.
//
// Part 3, colour (ColorFull): the soft-gluon emission density off a colour
//   singlet of partons is  -sum_{i<j} <T_i.T_j> W_ij  (on the sphere) or the
//   same with 2 K_ij (on the plane). ColorFull gives the colour correlators
//   <T_i.T_j> = <M|T_i.T_j|M> / <M|M>. For a singlet, colour conservation
//   sum_i T_i = 0 gives T_i.T_j = (C_k - C_i - C_j)/2 for three partons and
//   T_q.T_qbar = -CF for a q qbar dipole. We check those, then evaluate the
//   emission density for q qbar g on both sides of the projection.

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include "ColorFull/Core/Col_functions.h"
#include "ColorFull/Core/Col_amp.h"

using namespace ColorFull;
using namespace std;

typedef array<double, 2> Vec2; // point on the transverse plane
typedef array<double, 3> Vec3; // point on the sphere (or in R^3)

double dot(const Vec2 &a, const Vec2 &b) { return a[0] * b[0] + a[1] * b[1]; }
double dot(const Vec3 &a, const Vec3 &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec2 minus2(const Vec2 &a, const Vec2 &b) { return {a[0] - b[0], a[1] - b[1]}; }

// ----- the maps -------------------------------------------------------------

// plane -> unit sphere (notebook: rr[x] with r = 1)
Vec3 to_sphere(const Vec2 &x)
{
	double xx = dot(x, x);
	return {2 * x[0] / (1 + xx), 2 * x[1] / (1 + xx), (1 - xx) / (1 + xx)};
}

// unit sphere -> plane (inverse of to_sphere)
Vec2 to_plane(const Vec3 &n) { return {n[0] / (1 + n[2]), n[1] / (1 + n[2])}; }

// plane at height d -> sphere of radius r, projecting from s = (0,0,-r)
// (notebook: rrd[x])
Vec3 to_sphere_rd(const Vec2 &x, double r, double d)
{
	double xx = dot(x, x), D = d + r;
	return {2 * r * D * x[0] / (D * D + xx), 2 * r * D * x[1] / (D * D + xx),
		r * (D * D - xx) / (D * D + xx)};
}

// conformal factor: dOmega = Omega(x)^2 d^2x
double Omega(const Vec2 &x) { return 2 / (1 + dot(x, x)); }

// ----- kernels --------------------------------------------------------------

double n_ab(const Vec3 &a, const Vec3 &b) { return 1 - dot(a, b); }

// BMS: W_ab(k) = n_ab / (n_ak n_kb)
double W_bms(const Vec3 &a, const Vec3 &b, const Vec3 &k)
{
	return n_ab(a, b) / (n_ab(a, k) * n_ab(k, b));
}

// JIMWLK dipole kernel: (x-y)^2 / ((x-z)^2 (z-y)^2)
double K_jimwlk(const Vec2 &x, const Vec2 &y, const Vec2 &z)
{
	Vec2 xy = minus2(x, y), xz = minus2(x, z), zy = minus2(z, y);
	return dot(xy, xy) / (dot(xz, xz) * dot(zy, zy));
}

// <T_i.T_j> = <M|T_i.T_j|M> / <M|M> at the Nc currently set in cf
double TiTj(const Col_functions &cf, const Col_amp &M, int i, int j)
{
	return cf.double_num(cf.color_correlator(M, i, j)) /
	       cf.double_num(cf.scalar_product(M, M));
}

int main()
{
	mt19937 rng(2023);
	uniform_real_distribution<double> U(-3.0, 3.0);
	auto random_point = [&]()
	{ return Vec2{U(rng), U(rng)}; };
	const int N = 10000;

	cout << scientific << setprecision(2);

	// ======================= Part 1: the map =======================
	double err_norm = 0, err_inv = 0, err_tan = 0, err_rd = 0, err_line = 0;
	double err_conf = 0, err_ident = 0;
	const double r = 1.7, d = 0.4; // any sphere radius / plane height
	for (int i = 0; i < N; ++i)
	{
		Vec2 x = random_point(), y = random_point();
		Vec3 n = to_sphere(x), m = to_sphere(y);

		// |n| = 1
		err_norm = max(err_norm, fabs(dot(n, n) - 1));
		// to_plane(to_sphere(x)) = x
		Vec2 back = minus2(to_plane(n), x);
		err_inv = max(err_inv, sqrt(dot(back, back)));
		// |x| = tan(theta/2), theta = polar angle of n
		double theta = acos(n[2]);
		err_tan = max(err_tan, fabs(sqrt(dot(x, x)) - tan(theta / 2)));

		// general r, d: on the sphere of radius r, and on the line s + t (X - s)
		Vec3 w = to_sphere_rd(x, r, d);
		err_rd = max(err_rd, fabs(sqrt(dot(w, w)) - r));
		// the line from s=(0,0,-r) to X=(x1,x2,d): t fixed by the z-component
		double t = (w[2] + r) / (d + r);
		err_line = max(err_line, fabs(w[0] - t * x[0]) + fabs(w[1] - t * x[1]));

		// conformality: columns of the Jacobian dn/dx1, dn/dx2 (finite
		// differences) are orthogonal with length Omega
		const double h = 1e-6;
		Vec3 J[2];
		for (int c = 0; c < 2; ++c)
		{
			Vec2 xp = x, xm = x;
			xp[c] += h;
			xm[c] -= h;
			Vec3 np = to_sphere(xp), nm = to_sphere(xm);
			for (int k = 0; k < 3; ++k)
				J[c][k] = (np[k] - nm[k]) / (2 * h);
		}
		double O2 = Omega(x) * Omega(x);
		err_conf = max(err_conf, (fabs(dot(J[0], J[0]) - O2) + fabs(dot(J[1], J[1]) - O2) +
					  fabs(dot(J[0], J[1]))) /
					     O2);

		// 1 - n.m = 2 (x-y)^2 / ((1+x^2)(1+y^2))
		Vec2 xy = minus2(x, y);
		double rhs = 2 * dot(xy, xy) / ((1 + dot(x, x)) * (1 + dot(y, y)));
		err_ident = max(err_ident, fabs(n_ab(n, m) - rhs) / rhs);
	}
	cout << "Part 1: stereographic map, worst error over " << N << " random points\n"
	     << "  |n(x)| = 1                                  " << err_norm << "\n"
	     << "  to_plane(to_sphere(x)) = x                  " << err_inv << "\n"
	     << "  |x| = tan(theta/2)                          " << err_tan << "\n"
	     << "  rrd: |w| = r  (r = " << fixed << setprecision(1) << r << ", d = " << d
	     << ")" << scientific << setprecision(2) << "            " << err_rd << "\n"
	     << "  rrd: w on line from (0,0,-r) to (x,d)       " << err_line << "\n"
	     << "  conformal: J^T J = Omega^2 1  (relative)    " << err_conf
	     << "   (finite differences)\n"
	     << "  1-n.m = 2(x-y)^2/((1+x^2)(1+y^2)) (rel.)    " << err_ident << "\n\n";

	// ======================= Part 2: kernels =======================
	double rmin = 1e9, rmax = -1e9;
	for (int i = 0; i < N; ++i)
	{
		Vec2 xa = random_point(), xb = random_point(), z = random_point();
		Vec3 na = to_sphere(xa), nb = to_sphere(xb), nk = to_sphere(z);
		double sphere = W_bms(na, nb, nk) * Omega(z) * Omega(z); // W dOmega / d^2z
		double plane = K_jimwlk(xa, xb, z);
		double ratio = sphere / plane;
		rmin = min(rmin, ratio);
		rmax = max(rmax, ratio);
	}
	cout << fixed << setprecision(12)
	     << "Part 2: (W_ab dOmega_k) / (K_ab d^2z) over " << N << " random a, b, k\n"
	     << "  min = " << rmin << "\n  max = " << rmax << "   (expect exactly 2)\n\n";

	// ======================= Part 3: colour =======================
	Col_functions cf; // TR = 1/2 by default
	// q qbar singlet: quark 1, antiquark 2
	Col_amp dipole("[{1,2}]");
	// q g qbar singlet: quark 1, gluon 2, antiquark 3, i.e. (t^2)_{1 3}
	Col_amp qgqb("[{1,2,3}]");

	cout << defaultfloat << "Part 3: colour correlators from ColorFull\n"
	     << "  exact <q qbar|T_1.T_2|q qbar> = " << cf.color_correlator(dipole, 1, 2)
	     << "   (divide by <M|M> = " << cf.scalar_product(dipole, dipole) << ")\n";

	cout << fixed << setprecision(6);
	for (double Nc : {3.0, 5.0})
	{
		cf.set_Nc(Nc);
		double TR = cf.get_TR();
		double CA = 2 * TR * Nc, CF = TR * (Nc * Nc - 1) / Nc;

		double dip = TiTj(cf, dipole, 1, 2);
		double T_q_g = TiTj(cf, qgqb, 1, 2);
		double T_g_qb = TiTj(cf, qgqb, 2, 3);
		double T_q_qb = TiTj(cf, qgqb, 1, 3);

		cout << "  Nc = " << setprecision(0) << Nc << setprecision(6) << "  (CA = " << CA
		     << ", CF = " << CF << ")\n"
		     << "    q qbar:    T_q.T_qb = " << setw(10) << dip << "   expect -CF      = " << -CF << "\n"
		     << "    q g qbar:  T_q.T_g  = " << setw(10) << T_q_g << "   expect -CA/2    = " << -CA / 2 << "\n"
		     << "               T_g.T_qb = " << setw(10) << T_g_qb << "   expect -CA/2    = " << -CA / 2 << "\n"
		     << "               T_q.T_qb = " << setw(10) << T_q_qb << "   expect CA/2-CF  = " << CA / 2 - CF << "\n"
		     << "    sum_j T_q.T_j + T_q^2 = " << T_q_g + T_q_qb + CF << "   (colour conservation: 0)\n";

		// Emission density of a soft gluon at k off the q g qbar system.
		// Directions of q, g, qbar and the soft gluon, chosen on the plane:
		Vec2 xq{0.3, -0.2}, xg{-0.8, 0.5}, xqb{1.1, 0.9}, z{0.2, 0.7};
		Vec3 nq = to_sphere(xq), ng = to_sphere(xg), nqb = to_sphere(xqb), nk = to_sphere(z);

		double on_sphere = -(T_q_g * W_bms(nq, ng, nk) + T_g_qb * W_bms(ng, nqb, nk) +
				     T_q_qb * W_bms(nq, nqb, nk)) *
				   Omega(z) * Omega(z);
		double on_plane = -2 * (T_q_g * K_jimwlk(xq, xg, z) + T_g_qb * K_jimwlk(xg, xqb, z) +
					T_q_qb * K_jimwlk(xq, xqb, z));
		cout << "    emission density per d^2z, sphere side: " << on_sphere
		     << "\n                               plane side:  " << on_plane << "\n";
	}
	return 0;
}
