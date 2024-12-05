/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

//redox; based off of file /opt/openfoam9/src/thermophysicalModels/specie/reaction/Reactions/ReversibleReaction/ReversibleReaction.*

\*---------------------------------------------------------------------------*/

#include "ReversibleSolidHeterogeneousReaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionRate>
Foam::ReversibleSolidHeterogeneousReaction<ReactionRate>::ReversibleSolidHeterogeneousReaction
(
    const solidHeterogeneousReaction& reaction,
    const ReactionRate& k,
    const scalar nReact
)
:
    solidHeterogeneousReaction(reaction),
    k_(k),
	Kcs_(List<scalar>()), //redox
    heatReact_(List<scalar>()),
    nReact_(nReact)
{}


template<class ReactionRate>
Foam::ReversibleSolidHeterogeneousReaction<ReactionRate>::ReversibleSolidHeterogeneousReaction
(
    const PtrList<volScalarField>& gasPhaseGases,
    const speciesTable& components,
    Istream& is,
    const speciesTable& pyrolysisGases
)
:
    solidHeterogeneousReaction(gasPhaseGases, components, is, pyrolysisGases),
    k_(components, is),
	Kcs_(readScalar(is)), //redox
    heatReact_(readScalar(is)),
    nReact_(List<scalar>())
{
    List<int> lhs;
    DynamicList<scalar> dnReact;
    lhs.append(slhs());
    lhs.append(glhs());
    forAll(lhs,i)
    {
        dnReact.append(readScalar(is));
    }
	//redox
	List<int> rhs;
	rhs.append(srhs());
	rhs.append(grhs());
	forAll(rhs,i)
	{
		dnReact.append(readScalar(is));
	}
    is.readEnd("solidArrheniusReactionRate(Istream&)");
    nReact_ = dnReact.shrink();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionRate>
Foam::scalar Foam::ReversibleSolidHeterogeneousReaction<ReactionRate>::kf
(
    const scalar T,
    const scalar p,
	const scalar nonstoichiometry, //redox-delta
    const scalarField& c
) const
{  
	scalar kfactor;
	scalar prefactor;
	scalar x0, x1, f, df, aCoeff, bCoeff, cCoeff, dCoeff, xh2agent, xcoagent, xco, xco2;
	scalar Kmain, Kcomp, Xratio, Kceox;
	scalar goCehy, goCeco, Kcehy, Kceco;
	scalar xch4main, xch4comp, highDeltaStabilizer, xh2agentfactor, xcoagentfactor;
	
	x0 = 0.001;
	
	scalar goMain, goComp, goCeox;
	goMain = -196.65*T - 2.2808e4;
	goComp = 1.7438*T - 8.0210e5;
	goCehy = 56.0182*T - 2.4863e5;
	goCeco = 86.3547*T - 2.8204e5;
	goCeox = -86.2658*T + 2.8190e5;
	
	scalar pMultiple; //redox-pseudoKeq
	pMultiple = p/101325; //redox-pseudoKeq
	
	//for agent reaction rate control
	Kcehy = exp(-(goCehy+1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry))))/8.314/T);
	Kceco = exp(-(goCeco+1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry))))/8.314/T);

	//Newton-Raphson method - version 2
	
	//org
	//Kmain = exp(-(goMain+1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry))))/8.314/T);
	//Kcomp = exp(-(goComp+4*(1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry)))))/8.314/T);
	
	//with pseudoKeq
	Kmain = 1.0/pMultiple/pMultiple*exp(-(goMain+1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry))))/8.314/T);
	Kcomp = 1.0/pMultiple/pMultiple*exp(-(goComp+4*(1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry)))))/8.314/T);
	
	aCoeff = Kcomp/Kmain + 1.0;
	bCoeff = -3;
	cCoeff = 3;
	dCoeff = -1;
	f = aCoeff*pow(x0, 3) + bCoeff*pow(x0, 2) + cCoeff*x0 + dCoeff;
	df = 3*aCoeff*pow(x0, 2) + 2*bCoeff*x0 + cCoeff;
	x1 = x0 - f/df;
	// maximum 10 iterations
	for (int j = 1; j <= 10; j++) 
	{ 
		x0 = x1;
		f = aCoeff*pow(x0, 3) + bCoeff*pow(x0, 2) + cCoeff*x0 + dCoeff;
		df = 3*aCoeff*pow(x0, 2) + 2*bCoeff*x0 + cCoeff;
		x1 = x0 - f/df;
	}
	
	//XCO2 = (1-x1)/3;
	Xratio = (1-x1)/x1;

	//for agent stabilizers
	xh2agent = 1./(1.+Kcehy); //mol fraction of H2 for H2 agent reaction
	xcoagent = 1./(1.+Kceco); //mol fraction of CO for Co agent reaction

	//for oxidation
	Kceox = exp(-(goCeox-(1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry)))))/8.314/T);
	xco = Kceox/(1.+Kceox);
	xco2 = 1-xco;

	//for main stabilizer addition for high-delta low conversion
	aCoeff = 4.0;
	bCoeff = 0.;
	cCoeff = 3.0*Kmain;
	dCoeff = -1.*Kmain;
	f = aCoeff*pow(x0, 3) + bCoeff*pow(x0, 2) + cCoeff*x0 + dCoeff;
	df = 3*aCoeff*pow(x0, 2) + 2*bCoeff*x0 + cCoeff;
	x1 = x0 - f/df;
	// maximum 10 iterations
	for (int j = 1; j <= 10; j++) 
	{ 
		x0 = x1;
		f = aCoeff*pow(x0, 3) + bCoeff*pow(x0, 2) + cCoeff*x0 + dCoeff;
		df = 3*aCoeff*pow(x0, 2) + 2*bCoeff*x0 + cCoeff;
		x1 = x0 - f/df;
	}
	xch4main = 1.-3.*x1;//mol fraction of CH4 for main reaction
	
	//for comp stabilizer addition for high-delta low conversion
	//mol fraction of CH4 for comp reaction
	aCoeff = 4.0;
	bCoeff = 0.;
	cCoeff = 3.0*Kcomp;
	dCoeff = -1.*Kcomp;
	f = aCoeff*pow(x0, 3) + bCoeff*pow(x0, 2) + cCoeff*x0 + dCoeff;
	df = 3*aCoeff*pow(x0, 2) + 2*bCoeff*x0 + cCoeff;
	x1 = x0 - f/df;
	// maximum 10 iterations
	for (int j = 1; j <= 10; j++) 
	{ 
		x0 = x1;
		f = aCoeff*pow(x0, 3) + bCoeff*pow(x0, 2) + cCoeff*x0 + dCoeff;
		df = 3*aCoeff*pow(x0, 2) + 2*bCoeff*x0 + cCoeff;
		x1 = x0 - f/df;
	}
	xch4comp = 1.-3.*x1;//mol fraction of CH4 for comp reaction

	if (Kcs_ < 1.5) //syngas reaction, main reaction, wanted reaction
	{
		Xratio = Xratio * Xratio/5.;

		kfactor = 1/Xratio*(1/(1+1/Xratio));
		//kfactor = 1.;
		prefactor = 1.0;

		//check1
		if (kfactor>1.0)
		{
			kfactor = 1.0;
		}
		kfactor*=prefactor;
		
		//model validation
		highDeltaStabilizer = pow((1.-xch4main),0.5); //temporarily turned off
		//high-delta stabilizer 100%
		//highDeltaStabilizer = (1.-xch4main);
		
		kfactor*=highDeltaStabilizer;
		
		if (nonstoichiometry > 0.344)
		{	
			kfactor = 0.0;
		}
	}
	else if (Kcs_ < 2.5) //complete oxidation of CH4 reaction, initial reaction
	{
		Xratio = Xratio * Xratio/5.;

		//newton-raphson method-GOOD
		kfactor = 4*(1/(1+1/Xratio));
		prefactor = 1.0;
		kfactor*=prefactor;


		//v6-correct model validation
		highDeltaStabilizer = pow((1.-xch4comp),0.5); //temporarily turned off
		//high-delta stabilizer 100%
		//highDeltaStabilizer = (1.-xch4comp);

		kfactor*=highDeltaStabilizer;
		
		if (nonstoichiometry > 0.344)
		{	
			kfactor = 0.0;
		}
	}
	else if (Kcs_ < 3.5) //h2 agent, H2 + ceria = H2O + ceria
	{
		kfactor = 1.0;	
		prefactor = 1.0;
		kfactor*=prefactor;
	
		xh2agentfactor = xh2agent/(1.-xh2agent);

		kfactor =(0.8+2.2*1/(xh2agentfactor*xh2agentfactor/0.0015+1));

		//high-delta stabilizer-correct model validation
		kfactor *= pow((1.-xh2agent*xh2agent),0.5);
		//high-delta stabilizer 100p%
		//kfactor *= (1.-xh2agent);//kfactor *= (1.-xh2agent);

	}
	else if (Kcs_ < 4.5)//co agent, CO + ceria = CO2 + ceria
	{
		kfactor = 1.0;
		prefactor = 1.0;
		kfactor*=prefactor;

		xcoagentfactor =xcoagent/(1.-xcoagent);

		kfactor =(0.8+2.2*1/(xcoagentfactor*xcoagentfactor/0.0015+1));

		//high-delta stabilizer-correct model validation
		kfactor *= pow((1.-xcoagent*xcoagent),0.5);
		//high-delta stabilizer 100p%
		//kfactor *= (1.-xcoagent);//kfactor *= (1.-xcoagent);

	}
	else
	{

		kfactor = pow((1.-xco2),0.5);

		if (nonstoichiometry <0.001)
		{	
			kfactor = 0.0;
		}
	}
	
	//check2
	if (kfactor < 0.0)
	{
		kfactor = 3.5e-5;
	}
	//kfactor = 1.0; //uncomment if you don't want the stabilizer
	return k_(T, p, c)*kfactor;
}

template<class ReactionRate>
Foam::scalar Foam::ReversibleSolidHeterogeneousReaction<ReactionRate>::kr
(
	const scalar kfwd,
	const scalar T,
	const scalar p,
	const scalar nonstoichiometry, //redox-delta
	const scalarField & c
) const
{

	scalar krMax;
	krMax = 0.005;
	if (Kcs_ < 1.5) 
	{
		krMax = 0.005;
	}
	else if (Kcs_ < 2.5)
	{
		krMax = 0.005;
	}
	else if (Kcs_ < 3.5) 
	{
		krMax = 0.5;
	}
	else if (Kcs_ < 4.5)
	{
		krMax = 0.5;
	}
	else
	{
		krMax = 0.5;
	}

	//return Foam::min(kfwd/(this->Kcs(nonstoichiometry,T)),krMax);// krMax mode
	return kfwd/(this->Kcs(nonstoichiometry,T,p));// redox-pseduoKeq - added the p variable
}
	

//redox
template<class ReactionRate>
Foam::scalar Foam::ReversibleSolidHeterogeneousReaction<ReactionRate>::Kcs
(
	const scalar nonstoichiometry,
	const scalar T,
	const scalar p //redox-pseudoKeq
) const
{
	scalar KcsValue;
	scalar goMain, goComp, goCehy, goCeco, goCeox;
	scalar Kmain, Kcomp, MmixMain, MmixComp;
	scalar x0main, x1main, fMain, dfMain, aCoeff, bCoeff, cCoeff, dCoeff;
	scalar x0comp, x1comp, fComp, dfComp;
	scalar pMultiple; //redox-pseudoKeq

	goMain = -196.65*T - 2.2808e4;
	goComp = 1.7438*T - 8.0210e5;
	goCehy = 56.0182*T - 2.4863e5;
	goCeco = 86.3547*T - 2.8204e5;
	goCeox = -86.2658*T + 2.8190e5;

	//Newton-Raphson method - version 2
	x0main = 0.001;
	x0comp = 0.001;
	pMultiple = p/101325; //redox-pseudoKeq
	
	//org
	//Kmain = exp(-(goMain+1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry))))/8.314/T);
	//Kcomp = exp(-(goComp+4*(1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry)))))/8.314/T);
	
	//with pseudoKeq
	Kmain = 1.0/pMultiple/pMultiple*exp(-(goMain+1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry))))/8.314/T);
	Kcomp = 1.0/pMultiple/pMultiple*exp(-(goComp+4*(1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry)))))/8.314/T);	
	
	aCoeff = 4.;
	bCoeff = 0.;
	cCoeff = 3.;
	dCoeff = -1.;
	fMain = aCoeff/Kmain*pow(x0main, 3) + bCoeff*pow(x0main, 2) + cCoeff*x0main + dCoeff;
	fComp = aCoeff/Kcomp*pow(x0comp, 3) + bCoeff*pow(x0comp, 2) + cCoeff*x0comp + dCoeff;
	dfMain = 3*aCoeff/Kmain*pow(x0main, 2) + 2*bCoeff*x0main + cCoeff;
	dfComp = 3*aCoeff/Kcomp*pow(x0comp, 2) + 2*bCoeff*x0comp + cCoeff;
	x1main = x0main - fMain/dfMain;
	x1comp = x0comp - fComp/dfComp;
	// maximum 10 iterations
	for (int j = 1; j <= 10; j++) 
	{ 
		x0main = x1main;
		x0comp = x1comp;
		fMain = aCoeff/Kmain*pow(x0main, 3) + bCoeff*pow(x0main, 2) + cCoeff*x0main + dCoeff;
		fComp = aCoeff/Kcomp*pow(x0comp, 3) + bCoeff*pow(x0comp, 2) + cCoeff*x0comp + dCoeff;
		dfMain = 3*aCoeff/Kmain*pow(x0main, 2) + 2*bCoeff*x0main + cCoeff;
		dfComp = 3*aCoeff/Kcomp*pow(x0comp, 2) + 2*bCoeff*x0comp + cCoeff;
		x1main = x0main - fMain/dfMain;
		x1comp = x0comp - fComp/dfComp;
	}
	MmixMain = x1main*28. + 2.*x1main*2. + (1.-3*x1main)*16.;
	MmixComp = x1comp*44. + 2.*x1comp*18. + (1.-3*x1comp)*16.;
	
	if (Kcs_ < 1.5) //syngas reaction, main reaction, wanted reaction
	{	
		//org
		//KcsValue = 7.0/MmixMain/MmixMain*exp(-(goMain+1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry))))/8.314/T);
		
		//with pseudoKeq
		KcsValue = 1.0/pMultiple/pMultiple*7.0/MmixMain/MmixMain*exp(-(goMain+1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry))))/8.314/T);
		
	}

	else if (Kcs_ < 2.5) //complete oxidation of CH4 reaction, initial reaction
	{
		//org
		//KcsValue = 891.0/MmixComp/MmixComp*exp(-(goComp+4*(1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry)))))/8.314/T);
		
		//pseudoKeq
		KcsValue = 1.0/pMultiple/pMultiple*891.0/MmixComp/MmixComp*exp(-(goComp+4*(1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry)))))/8.314/T);
			
	}
	else if (Kcs_ < 3.5) //h2 agent, H2 + ceria = H2O + ceria
	{
		KcsValue = 9.0*exp(-(goCehy+1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry))))/8.314/T);

	}
	else if (Kcs_ < 4.5)//co agent, CO + ceria = CO2 + ceria
	{
		KcsValue = 1.57*exp(-(goCeco+1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry))))/8.314/T);

	}
	else
	{
		KcsValue = 0.636*exp(-(goCeox-(1000*(395.0-31.4*log10(nonstoichiometry))-T*(160+2.9*8.314*(log(0.345-nonstoichiometry)-log(nonstoichiometry)))))/8.314/T);
	}
	
	return KcsValue;
}


template<class ReactionRate>
Foam::scalar Foam::ReversibleSolidHeterogeneousReaction<ReactionRate>::heatReact() const
{
    return heatReact_;
}

template<class ReactionRate>
Foam::List<scalar> Foam::ReversibleSolidHeterogeneousReaction<ReactionRate>::nReact() const
{
	return nReact_;
}

template<class ReactionRate>
void Foam::ReversibleSolidHeterogeneousReaction<ReactionRate>::write
(
    Ostream& os
) const
{
    solidHeterogeneousReaction::write(os);
    os  << token::SPACE << "Reaction order: " << nReact_ << nl 
        << token::SPACE << "Energy of reaction (if applicable): " << heatReact_ << " J/kg" 
        << token::SPACE << nl << k_;
}


// ************************************************************************* //
