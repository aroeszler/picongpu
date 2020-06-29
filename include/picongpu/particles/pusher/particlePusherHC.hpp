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

/* used references:
 *
 * Introduction of the particle pusher: https://dx.doi.org/10.1063/1.4979989 or https://arxiv.org/abs/1701.05605 https://doi.org/10.3847/1538-4365/aab114
 *
 * Including little correction(s) by WarpX team: https://github.com/ECP-WarpX/WarpX/issues/320
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
        float_X const mass = attribute::getMass( weighting , particle );
        float_X const charge = attribute::getCharge( weighting , particle );

        using MomType = momentum::type;
        MomType const mom = particle[ momentum_ ];

        auto bField  = functorBField(pos);
        auto eField  = functorEField(pos);

        const float_X deltaT = DELTA_T;


        Gamma gamma;        
        
        const MomType mom_minus = mom + float_X(0.5) * charge * eField * deltaT;

        const sqrt_HC::float_X gamma_minus = gamma( mom_minus , mass );
        
        const sqrt_HC::float3_X tau = precisionCast<sqrt_HC::float_X>( float_X(0.5) * bField * charge * deltaT );

        const sqrt_HC::float_X sigma = pmacc::math::abs2( gamma_minus ) - pmacc::math::abs2( tau );
        
        const sqrt_HC::float_X u_star = pmacc::math::dot( precisionCast<sqrt_HC::float_X>( mom_minus ), tau ) / precisionCast<sqrt_HC::float_X>( mass * SPEED_OF_LIGHT );
       
        const sqrt_HC::float_X gamma_plus = math::sqrt( ( sigma + math::sqrt( pmacc::math::abs2( sigma ) + ( sqrt_HC::float_X(4.0) ) * ( pmacc::math::abs2( tau ) + pmacc::math::abs2( u_star ) ) ) ) * ( sqrt_HC::float_X(0.5) ) );
        
        const sqrt_HC::float3_X t_vector = ( ( sqrt_HC::float_X(1.0) ) / gamma_plus ) * tau;
        
        const sqrt_HC::float_X s = ( sqrt_HC::float_X(1.0) ) / ( sqrt_HC::float_X(1.0) + pmacc::math::abs2( t_vector ) );
                                    
        const MomType mom_plus = precisionCast<float_X>( s * ( precisionCast<sqrt_HC::float_X>( mom_minus ) + t_vector * pmacc::math::dot( precisionCast<sqrt_HC::float_X>( mom_minus ) , t_vector ) ) + ( pmacc::math::cross( precisionCast<sqrt_HC::float_X>( mom_minus ) , t_vector)));
                                                        
        const MomType new_mom = mom_plus + eField * charge * deltaT * float_X(0.5) + pmacc::math::cross( mom_plus , precisionCast<float_X>( t_vector ) );
         
        const sqrt_HC::float_X gamma_final = gamma( new_mom , mass );
                     
        particle[ momentum_ ] = new_mom;

        Velocity velocity;
        
        const float3_X vel = velocity( new_mom , mass );

        for( uint32_t d=0 ; d<simDim ; ++d )
            {
               pos[d] += ( vel[d] * deltaT ) / cellSize[d];
            }
        
        
    }

    static pmacc::traits::StringProperty getStringProperties()
    {
        pmacc::traits::StringProperty propList( "name", "HC" );
        return propList;
    }
};
} // namespace particlePusherHC
} // namespace picongpu
