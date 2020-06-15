/* Copyright 2013-2020 Heiko Burau, Rene Widera, Richard Pausch, Annegret Roeszler
 *
 * This file is part of PIConGPU.
 *
 * PIConGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PIConGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PIConGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "picongpu/simulation_defines.hpp"
#include "picongpu/traits/attribute/GetMass.hpp"
#include "picongpu/traits/attribute/GetCharge.hpp"

namespace picongpu
{
namespace particlePusherHC
{

template<class Velocity, class Gamma>
struct Push
{
    /* this is an optional extension for sub-sampling pushes that enables grid to particle interpolation
     * for particle positions outside the super cell in one push
     */
    using LowerMargin = typename pmacc::math::CT::make_Int<simDim,0>::type;
    using UpperMargin = typename pmacc::math::CT::make_Int<simDim,0>::type;

    template< typename T_FunctorFieldE, typename T_FunctorFieldB, typename T_Particle, typename T_Pos >
    HDINLINE void operator()(
        const T_FunctorFieldB functorBField,
        const T_FunctorFieldE functorEField,
        T_Particle & particle,
        T_Pos & pos,
        const uint32_t
    )
    {
        float_X const weighting = particle[ weighting_ ];
        float_X const mass = attribute::getMass( weighting, particle );
        float_X const charge = attribute::getCharge( weighting, particle );

        using MomType = momentum::type;
        MomType const mom = particle[ momentum_ ];

        auto bField  = functorBField(pos);
        auto eField  = functorEField(pos);

        const float_X QoM = charge / mass;

        const float_X deltaT = DELTA_T;


        Gamma gamma;
        
        /* ORIGINAL BORIS IMPLEMENTATION
        const float_X gamma_reci = float_X(1.0) / gamma(mom_minus, mass);
        const float3_X t = float_X(0.5) * QoM * bField * gamma_reci * deltaT;
        auto s  = float_X(2.0) * t * (float_X(1.0) / (float_X(1.0) + pmacc::math::abs2(t)));

        const MomType mom_prime = mom_minus + pmacc::math::cross(mom_minus, t);
        const MomType mom_plus = mom_minus + pmacc::math::cross(mom_prime, s);

        const MomType new_mom = mom_plus + float_X(0.5) * charge * eField * deltaT;

        particle[ momentum_ ] = new_mom;

        Velocity velocity;
        const float3_X vel = velocity(new_mom, mass);

        for(uint32_t d=0;d<simDim;++d)
        {
            pos[d] += (vel[d] * deltaT) / cellSize[d];
        }
        */
        
        // gleich überall in mom umrechnen
        //wie heißt das "alte" Momentum?
        
        //const sqrt_HC::float_X gamma_i = gamma( mom , mass ); //brauchen wir das überhaupt?
        
        //const sqrt_HC::float3_X x_first_half = pos[d] + mom * (deltaT * sqrt_HC::float_X(0.5)/(gamma_i * mass));
        
        const MomType mom_minus = mom + sqrt_HC::float_X(0.5) * charge * eField * deltaT;
        
        const sqrt_HC::float_X gamma_minus = gamma( mom_minus , mass );
        
        const sqrt_HC::float3_X tau = (sqrt_HC::float_X(0.5)) * bField * charge * deltaT;
        
        const sqrt_HC::float_X sigma = pmacc::math::abs2(gamma_minus) - pmacc::math::abs(tau); //ist die abs2 Fkt sowohl für float_X als auch float3_X passend?
        
        const sqrt_HC::float3_X u_star = pmacc::math::dot( (mom_minus * sqrt_HC::float_X(1.0) / mass ) , tau * ((sqrt_HC::float_X(1.0)) / c ));
        // müssen die alle const sein?
        
        const sqrt_HC::float_X gamma_plus = math::sqrt(( sigma + math::sqrt( pmacc::math::abs2(sigma) + (sqrt_HC::float_X(4.0)) * (pmacc::math::abs2(tau) + pmacc::math::abs2(u_star))))* (sqrt_HC::float_X(0.5)));
        
        const sqrt_HC::float3_X t = ((sqrt_HC::float_X(1.0))/gamma_plus) * tau;
        
        const sqrt_HC::float_X s = (sqrt_HC::float_X(1.0))/(sqrt_HC::float_X(1.0) + pmacc::math::abs2(t));
                                    
         const MomType mom_plus = (s * (mom_minus + t_vector * pmacc::math::dot( mom_minus , t_vector )) + (pmacc::math::cross(mom_minus ,t_vector))); //nochmal angucken, sollte man mom_minus/m nehmen und am ende *m nehmen? eher nicht
                                                        
         const MomType new_mom = mom_plus + eField * charge * deltaT * (sqrt_HC::float_X(0.5)) + pmacc::math::cross(mom_plus,t_vector);
        //geht das so einfach?
         
        const sqrt_HC::float_X gamma_final = gamma( new_mom, mass); //???
                     
         //position update
         //const sqrt_HC::float3_X x_i = x_first_half + new_mom * deltaT * ((sqrt_HC::float_X(0.5))/ (mass * gamma_final))
         particle[ momentum_ ] = new_mom;

         Velocity velocity;
        
         const float3_X vel = velocity(new_mom, mass);

         for(uint32_t d=0;d<simDim;++d)
            {
                pos[d] += (vel[d] * deltaT) / cellSize[d];
            }
        
        
    }

    /*static pmacc::traits::StringProperty getStringProperties()
    {
        pmacc::traits::StringProperty propList( "name", "Boris" );
        return propList;
    }*/
};
} // namespace particlePusherBoris
} // namespace picongpu
