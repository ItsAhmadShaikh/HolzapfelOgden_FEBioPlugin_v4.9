#include "stdafx.h"
#include "FEBioMech/FEElasticFiberMaterialUC.h"
#include "Holzapfel-Ogden.h"
#include <iostream>
#include <fstream>


//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(HolzapfelMyocardiumPI, FEUncoupledMaterial)
	ADD_PARAMETER(m_a, FE_RANGE_GREATER_OR_EQUAL(0.0), "a");
	ADD_PARAMETER(m_b, FE_RANGE_GREATER_OR_EQUAL(0.0), "b");
	ADD_PARAMETER(m_af, FE_RANGE_GREATER_OR_EQUAL(0.0), "af");
	ADD_PARAMETER(m_bf, FE_RANGE_GREATER_OR_EQUAL(0.0), "bf");
	ADD_PARAMETER(m_as, FE_RANGE_GREATER_OR_EQUAL(0.0), "as");
	ADD_PARAMETER(m_bs, FE_RANGE_GREATER_OR_EQUAL(0.0), "bs");
	ADD_PARAMETER(m_afs, FE_RANGE_GREATER_OR_EQUAL(0.0), "afs");
	ADD_PARAMETER(m_bfs, FE_RANGE_GREATER_OR_EQUAL(0.0), "bfs");
	ADD_PARAMETER(m_asn, FE_RANGE_GREATER_OR_EQUAL(0.0), "asn");
	ADD_PARAMETER(m_bsn, FE_RANGE_GREATER_OR_EQUAL(0.0), "bsn");
	ADD_PARAMETER(m_anf, FE_RANGE_GREATER_OR_EQUAL(0.0), "anf");
	ADD_PARAMETER(m_bnf, FE_RANGE_GREATER_OR_EQUAL(0.0), "bnf");
END_FECORE_CLASS();

mat3ds HolzapfelMyocardiumPI::sbar(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	//deformation gradient (deviatoric part)
	double     J = pt.m_J;
	mat3d   Fbar = pt.m_F* pow(J, -1.0 / 3.0);
	mat3ds  bbar = pt.DevLeftCauchyGreen();
	mat3ds  CRbar = pt.DevRightCauchyGreen();
	
	// get the initial fiber direction
    mat3d Q = GetLocalCS(mp);
    
	vec3d e10;
	e10.x = Q[0][0];
	e10.y = Q[1][0];
	e10.z = Q[2][0];

	vec3d e20;
	e20.x = Q[0][1];
	e20.y = Q[1][1];
	e20.z = Q[2][1];

	vec3d e30;
	e30.x = Q[0][2];
	e30.y = Q[1][2];
	e30.z = Q[2][2];

	//Spatial material axes
	vec3d e1 = Fbar*e10;
	vec3d e2 = Fbar*e20;
	vec3d e3 = Fbar*e30;

	
	// calculate the current FIBER, SHEET, NORMAL vectors
	vec3d ft = e1;	// current FIBER 
	vec3d st = e2;	// current SHEET 
	vec3d nt = e3;	// current NORMAL
	
	// calculate the outer product of ft & st
	mat3ds fxf = dyad(ft);
	mat3ds sxs = dyad(st);
	mat3d  fxs = ft & st;
	mat3d  sxn = st & nt;
	mat3d  nxf = nt & ft;

	mat3ds fxs_sym = fxs.sym();
	mat3ds sxn_sym = sxn.sym();
	mat3ds nxf_sym = nxf.sym();

	// Isochoric Invariants
	double I1   = bbar.tr();
	double I4f = e10*(CRbar*e10);
	double I4s = e20*(CRbar*e20);
	
	double I8fs = e10*(CRbar*e20);
	double I8sn = e20*(CRbar*e30);
	double I8nf = e30*(CRbar*e10);

	//consider fibers under tension, only
	if (I4f    < 1.0)
	{
		I4f = 1.0;
	}
	if (I4s    < 1.0)
	{
		I4s = 1.0;
	}
    
    // Get Material Properties
	double a = m_a(mp);
	double b = m_b(mp);
	double af = m_af(mp);
	double bf = m_bf(mp);
	double as = m_as(mp);
	double bs = m_bs(mp);
	double afs = m_afs(mp);
	double bfs = m_bfs(mp);
	double asn = m_asn(mp);
	double bsn = m_bsn(mp);
	double anf = m_anf(mp);
	double bnf = m_bnf(mp);
   
	//Strain Energy derivatives (isochoric part)
	double W1 = (a / 2.0)*exp(b*(I1 - 3.0));
	double W4f = (af*(I4f - 1.0))*exp(bf*pow(I4f - 1.0, 2.0));
	double W4s = (as*(I4s - 1.0))*exp(bs*pow(I4s - 1.0, 2.0)); 
	
	double W8fs = (afs*I8fs)*exp(bfs*pow(I8fs, 2.0)); 
	double W8sn = (asn*I8sn)*exp(bsn*pow(I8sn, 2.0));
	double W8nf = (anf*I8nf)*exp(bnf*pow(I8nf, 2.0));

	mat3ds tbar_coupling = ((2.0*W8fs)*(fxs_sym))+((2.0*W8sn)*(sxn_sym))+((2.0*W8nf)*(nxf_sym));
	
	mat3ds tbar = ((2.0*W1)*bbar) + ((2.0*W4f)*fxf) + ((2.0*W4s)*sxs) + (tbar_coupling); //Last term needs 2.0 factor because of .sym()=0.5(A+At)
	
	mat3ds sbar = pow(J, -1.0)*tbar;
	
	return sbar;
	
}

//-----------------------------------------------------------------------------------------------------------------------
tens4ds HolzapfelMyocardiumPI::cbar(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	//deformation gradient (deviatoric part)
	double     J = pt.m_J;
	mat3d   Fbar = pt.m_F* pow(J, -1.0 / 3.0);
	mat3ds  bbar = pt.DevLeftCauchyGreen();
	mat3ds  CRbar = pt.DevRightCauchyGreen();
	
	// get the initial fiber direction
    mat3d Q = GetLocalCS(mp);
    
	vec3d e10;
	e10.x = Q[0][0];
	e10.y = Q[1][0];
	e10.z = Q[2][0];

	vec3d e20;
	e20.x = Q[0][1];
	e20.y = Q[1][1];
	e20.z = Q[2][1];

	vec3d e30;
	e30.x = Q[0][2];
	e30.y = Q[1][2];
	e30.z = Q[2][2];

	//current material axes
	vec3d e1 = Fbar*e10;
	vec3d e2 = Fbar*e20;
	vec3d e3 = Fbar*e30;


	// calculate the current FIBER, SHEET, NORMAL vectors
	vec3d ft = e1;	// current FIBER 
	vec3d st = e2;	// current SHEET 
	vec3d nt = e3;	// current NORMAL

	// calculate the outer product of ft & st
	mat3ds fxf = dyad(ft);
	mat3ds sxs = dyad(st);
	mat3d  fxs = ft & st;
	mat3d  sxn = st & nt;
	mat3d  nxf = nt & ft;

	mat3ds fxs_sym = fxs.sym();
	mat3ds sxn_sym = sxn.sym();
	mat3ds nxf_sym = nxf.sym();

	//Invariants
	double I1  = bbar.tr();
	double I4f = e10*(CRbar*e10);
	double I4s = e20*(CRbar*e20);
	
	double I8fs = e10*(CRbar*e20);
	double I8sn = e20*(CRbar*e30);
	double I8nf = e30*(CRbar*e10);
	
	//consider fibers under tension, only
	if (I4f    < 1.0)
	{
		I4f = 1;
	}
	if (I4s    < 1.0)
	{
		I4s = 1;
	}
	
    // Get Material Properties
	double a = m_a(mp);
	double b = m_b(mp);
	double af = m_af(mp);
	double bf = m_bf(mp);
	double as = m_as(mp);
	double bs = m_bs(mp);
	double afs = m_afs(mp);
	double bfs = m_bfs(mp);
	double asn = m_asn(mp);
	double bsn = m_bsn(mp);
	double anf = m_anf(mp);
	double bnf = m_bnf(mp);    
    
	//Strain Energy second derivatives (isochoric part)
	double W1_dot = a*b / 2.0*exp(b*(I1 - 3));
	double W4f_dot = af*(2.0*bf*pow((I4f - 1.0), 2.0) + 1)*exp(bf*pow((I4f - 1), 2.0));
	double W4s_dot = as*(2.0*bs*pow((I4s - 1.0), 2.0) + 1)*exp(bs*pow((I4s - 1), 2.0));
	
	double W8fs_dot = afs*(2.0*bfs*pow(I8fs, 2.0) + 1)*exp(bfs*pow(I8fs, 2.0));
	double W8sn_dot = asn*(2.0*bsn*pow(I8sn, 2.0) + 1)*exp(bsn*pow(I8sn, 2.0));
	double W8nf_dot = anf*(2.0*bnf*pow(I8nf, 2.0) + 1)*exp(bnf*pow(I8nf, 2.0));
	
	//Calculate cbar
	tens4ds cbar_coupling = W8fs_dot*dyad1s(2.0*fxs_sym)+ W8sn_dot*dyad1s(2.0*sxn_sym)+ W8nf_dot*dyad1s(2.0*nxf_sym);
	
	tens4ds cbar = (4.0*W1_dot)*dyad1s(bbar) + (4.0*W4f_dot)*dyad1s(fxf) + (4.0*W4s_dot)*dyad1s(sxs) + cbar_coupling; 
	return cbar;
}

//----------------------------------------------------------------------------------------------------------------
mat3ds HolzapfelMyocardiumPI::DevStress(FEMaterialPoint& mp)
{
	mat3ds sbar = HolzapfelMyocardiumPI::sbar(mp);
	return sbar.dev();
}


//----------------------------------------------------------------------------------------------------------------
tens4ds HolzapfelMyocardiumPI::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3ds  I(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);

	double  J = pt.m_J;
	mat3ds  sbar = HolzapfelMyocardiumPI::sbar(mp);
	mat3ds  siso = sbar.dev();                      

	tens4ds IxI = dyad1s(I);
	tens4ds I4 = dyad4s(I);
	tens4ds sisoxI = dyad1s(siso, I);

	tens4ds proj = I4 - (1.0 / 3.0)*IxI;

	tens4ds cbarJi = (1.0 / J)*HolzapfelMyocardiumPI::cbar(mp);
	tens4ds proj_cbarJi_proj = cbarJi - (1.0 / 3.0)*ddots(cbarJi, IxI) + (1.0 / 9.0)*IxI*cbarJi.tr();

	return  proj_cbarJi_proj + (((2.0 / 3.0)*sbar.tr())*proj - (2.0 / 3.0)*sisoxI);
	
}


//-----------------------------------------------------------------------------------------------------------
double HolzapfelMyocardiumPI::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	//deformation gradient (deviatoric part)
	double     J = pt.m_J;
	mat3d   Fbar = pt.m_F*pow(J, -1.0 / 3.0);
	mat3ds  bbar = pt.DevLeftCauchyGreen();
	mat3ds  CRbar = pt.DevRightCauchyGreen();
	
	// get the initial fiber direction
    mat3d Q = GetLocalCS(mp);
    
	vec3d e10;
	e10.x = Q[0][0];
	e10.y = Q[1][0];
	e10.z = Q[2][0];

	vec3d e20;
	e20.x = Q[0][1];
	e20.y = Q[1][1];
	e20.z = Q[2][1];

	vec3d e30;
	e30.x = Q[0][2];
	e30.y = Q[1][2];
	e30.z = Q[2][2];

	//current material axes
	vec3d e1 = Fbar*e10;
	vec3d e2 = Fbar*e20;
	vec3d e3 = Fbar*e30;
	
	// calculate the current FIBER, SHEET, NORMAL vectors
	vec3d ft = e1;	// current FIBER 
	vec3d st = e2;	// current SHEET 
	vec3d nt = e3;	// current NORMAL

	// calculate the outer product of ft & st
	mat3ds fxf = dyad(ft);
	mat3ds sxs = dyad(st);
	mat3d  fxs = ft & st;
	mat3d  sxn = st & nt;
	mat3d  nxf = nt & ft;

	mat3ds fxs_sym = fxs.sym();
	mat3ds sxn_sym = sxn.sym();
	mat3ds nxf_sym = nxf.sym();
	

	//Invariants
	double I1  = bbar.tr();
	double I4f = e10*(CRbar*e10);
	double I4s = e20*(CRbar*e20);
	
	double I8fs = e10*(CRbar*e20);
	double I8sn = e20*(CRbar*e30);
	double I8nf = e30*(CRbar*e10);

	//consider fibers under tension, only
	if (I4f < 1.0)
	{
		I4f = 1.0;
	}
	if (I4s < 1.0)
	{
		I4s = 1.0;
	}	
    
    // Get Material Properties
	double a = m_a(mp);
	double b = m_b(mp);
	double af = m_af(mp);
	double bf = m_bf(mp);
	double as = m_as(mp);
	double bs = m_bs(mp);
	double afs = m_afs(mp);
	double bfs = m_bfs(mp);
	double asn = m_asn(mp);
	double bsn = m_bsn(mp);
	double anf = m_anf(mp);
	double bnf = m_bnf(mp);
	
	// Initialize
	double sedI1 = 0.0;
	double sedI4f = 0.0;
	double sedI4s = 0.0;
	double sedI8fs = 0.0;
	double sedI8sn = 0.0;
	double sedI8nf = 0.0;
	
	// Isotropic term
	if (b > 0.0)
	{
		sedI1 = (a / (2.0*b))*(exp(b*(I1 - 3.0))-1.0);
	
	} else {
	
		sedI1 = (a / 2.0)*(I1 - 3.0);
	
	}
	
	// Fiber term
	if (bf > 0.0)
	{
		
		sedI4f = (af / (2.0*bf))*(exp(bf*pow((I4f - 1.0), 2.0)) - 1.0);
		
	} else {
		
		sedI4f = (af / 2.0)*pow((I4f - 1.0), 2.0);
		
	}
		
	// Sheet term
	if (bs > 0.0)
	{
		
		sedI4s = (as / (2.0*bs))*(exp(bs*pow((I4s - 1.0), 2.0)) - 1.0);
		
	} else {
		
		sedI4s = (as / 2.0)*pow((I4s - 1.0), 2.0);
		
	}
	
	// FS coupling term
	if (bfs > 0.0)
	{
		
		sedI8fs = (afs / (2.0*bfs))*(exp(bfs*pow((I8fs), 2.0)) - 1.0);
		
	} else {
		
		sedI8fs = (afs / 2.0)*pow((I8fs), 2.0);
		
	}
	
	// SN coupling term
	if (bsn > 0.0)
	{
		
		sedI8sn = (asn / (2.0*bsn))*(exp(bsn*pow((I8sn), 2.0)) - 1.0);
		
	} else {
		
		sedI8sn = (asn / 2.0)*pow((I8sn), 2.0);
		
	}
	
	// NF coupling term
	if (bnf > 0.0)
	{
		
		sedI8nf = (anf / (2.0*bnf))*(exp(bnf*pow((I8nf), 2.0)) - 1.0);
		
	} else {
		
		sedI8nf = (anf / 2.0)*pow((I8nf), 2.0);
		
	}

	// Sum every component
	double sed = sedI1 + sedI4f + sedI4s + sedI8fs + sedI8sn + sedI8nf;

	return sed;
}


