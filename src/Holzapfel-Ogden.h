#pragma once
#include "stdafx.h"
#include "FEBioMech/FEUncoupledMaterial.h"
#include "FEBioMech/FEElasticFiberMaterialUC.h"

//-----------------------------------------------------------------------------
//! Holzapfel Myocardium Material

//! Implementation of a  Holzapfel - Ogden  uncoupled hyperelastic material.
class HolzapfelMyocardiumPI : public FEUncoupledMaterial
{
public:
	HolzapfelMyocardiumPI(FEModel* pfem) : FEUncoupledMaterial(pfem) {}

public:
	// a's : Stress Units
	// b's : Dimensionless parameters
 	FEParamDouble	m_a;			//! a :  Isotropic term parameter
	FEParamDouble	m_b;			//! b :  Isotropic term parameter
	FEParamDouble	m_af;			//! af:  Muscle Fiber tension parameter 
	FEParamDouble	m_bf;			//! bf:  Muscle Fiber tension parameter
	FEParamDouble	m_as;			//! as:  Collagen Fiber tension parameter 
	FEParamDouble	m_bs;			//! bs:  Collagen Fiber tension parameter
	FEParamDouble	m_afs;		//! afs:  Orthotropic term parameter
	FEParamDouble	m_bfs;		//! bfs:  Orthotropic term parameter
	FEParamDouble	m_asn;		//! Additional coupling parameter
	FEParamDouble	m_bsn;		//! Additional coupling parameter
	FEParamDouble	m_anf;		//! Additional coupling parameter
	FEParamDouble	m_bnf;		//! Additional coupling parameter
	

public:
	//! Deviatoric Cauchy stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& mp);

	//! Deviatoric spatial tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& mp);

	//! Devatioric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp);

	//! explicitly calculate the projection of tau_bar
	mat3ds sbar(FEMaterialPoint& mp);

	//! explicitly calculate the projection of c_bar
	tens4ds cbar(FEMaterialPoint& mp);

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
