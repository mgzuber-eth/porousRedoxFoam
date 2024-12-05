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

Class
    Foam::ReversibleSolidHeterogeneousReaction

Description
    Simple extension of Reaction to handle reversible reactions

SourceFiles
    ReversibleSolidHeterogeneousReaction.C

\*---------------------------------------------------------------------------*/

#ifndef ReversibleSolidHeterogeneousReaction_H
#define ReversibleSolidHeterogeneousReaction_H

#include "solidHeterogeneousReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ReversibleSolidHeterogeneousReaction Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionRate>
class ReversibleSolidHeterogeneousReaction
:
    public solidHeterogeneousReaction
{
    // Private data

        //- Reaction rate
        ReactionRate k_;

		//- Equilibrium constant
        scalar Kcs_; //redox
	
        //- Reaction heat
        scalar heatReact_;

        //- Reaction order
        List<scalar> nReact_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=
        (
            const ReversibleSolidHeterogeneousReaction<ReactionRate>&
        );


public:

    //- Runtime type information
    TypeName("reversible");


    // Constructors

        //- Construct from components
        ReversibleSolidHeterogeneousReaction
        (
            const solidHeterogeneousReaction& reaction,
            const ReactionRate& k,
            const scalar nReact
        );


        //- Construct from Istream
        ReversibleSolidHeterogeneousReaction
        (
	        const PtrList<volScalarField>& gasPhaseGases,
            const speciesTable& components,
            Istream& is,
            const speciesTable& pyrolysisGases
        );


    //- Destructor
    virtual ~ReversibleSolidHeterogeneousReaction()
    {}


    // Member Functions

        // ReversibleSolidHeterogeneousReaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar T,
                const scalar p,
				const scalar nonstoichiometry, //redox-delta
                const scalarField& c
            ) const;
	
			//- Reverse rate constant
			virtual scalar kr
			(
				const scalar kfwd,
				const scalar T,
				const scalar p,
				const scalar nonstoichiometry, //redox-delta
				const scalarField& c
			) const;

			//- Equilibrium constant - redox
            virtual scalar Kcs
			(
				const scalar nonstoichiometry, //redox-delta
				const scalar T,
				const scalar p //redox-pseudoKeq
			) const;

           //- Reaction order
            virtual scalar heatReact() const;

           //- Reaction order
           virtual List<scalar> nReact() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReversibleSolidHeterogeneousReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
