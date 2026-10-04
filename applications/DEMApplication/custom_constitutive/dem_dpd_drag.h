// Surface-based DPD--DEM constitutive law for suspended DEM particles.
//
// The DPD interaction coordinate is the signed distance from the DEM surface:
//
//     d_s = ||x_1 - x_2|| - R_DEM
//
// rather than the center-to-center distance. The particle identified as the
// DEM/suspended particle must have the larger interaction radius. In a DPD--DEM
// pair, this law uses:
//
//     d_s = r_center - max(R_1, R_2)
//
// so it is independent of the ordering of the two particles in the pair loop.
//
// The DPD force is active within:
//
//     d_s < DPD_CUTOFF_RADIUS.
//
// Positive d_s means the DPD-particle center lies outside the DEM surface;
// negative d_s means penetration into the DEM particle. Negative d_s is
// intentionally kept: it produces a stronger outward conservative force and
// expels penetrated DPD particles.

#if !defined(KRATOS_DEM_DPD_DRAG_H_INCLUDED)
#define KRATOS_DEM_DPD_DRAG_H_INCLUDED

#include "includes/define.h"
#include "includes/serializer.h"
#include "custom_constitutive/DEM_D_Linear_viscous_Coulomb_CL.h"

namespace Kratos
{

    class SphericParticle;

    /**
     * @class DEM_DPD_DRAG
     * @brief Surface-distance DPD--DEM law for a finite-radius suspended particle.
     *
     * Use this law only for DPD--suspended-DEM material pairs. Do not assign it to
     * DEM--DEM pairs, because DEM--DEM contact should retain its normal overlap-
     * based contact law.
     */
    class KRATOS_API(DEM_APPLICATION) DEM_DPD_DRAG
        : public DEM_D_Linear_viscous_Coulomb
    {
    public:
        using DEMDiscontinuumConstitutiveLaw::CalculateNormalForce;

        KRATOS_CLASS_POINTER_DEFINITION(DEM_DPD_DRAG);

        DEM_DPD_DRAG() = default;
        ~DEM_DPD_DRAG() override = default;

        std::string GetTypeOfLaw() override;

        /**
         * DPD coupling has a finite range independent of DEM overlap.
         */
        bool IsRangeForce() const override
        {
            return true;
        }

        void Check(Properties::Pointer pProp) const override;

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        std::unique_ptr<DEMDiscontinuumConstitutiveLaw> CloneUnique() override;

        void InitializeContact(
            SphericParticle *const element1,
            SphericParticle *const element2,
            const double indentation) override;

        void InitializeContactWithFEM(
            SphericParticle *const element,
            Condition *const wall,
            const double indentation,
            const double ini_delta = 0.0) override;

        void CalculateForces(
            const ProcessInfo &rProcessInfo,
            const double OldLocalElasticContactForce[3],
            double LocalElasticContactForce[3],
            double LocalDeltDisp[3],
            double LocalRelVel[3],
            double indentation,
            double previous_indentation,
            double ViscoDampingLocalContactForce[3],
            double &rCohesiveForce,
            SphericParticle *element1,
            SphericParticle *element2,
            bool &rSliding,
            double LocalCoordSystem[3][3]) override;

        void CalculateForcesWithFEM(
            const ProcessInfo &rProcessInfo,
            const double OldLocalElasticContactForce[3],
            double LocalElasticContactForce[3],
            double LocalDeltDisp[3],
            double LocalRelVel[3],
            double indentation,
            double previous_indentation,
            double ViscoDampingLocalContactForce[3],
            double &rCohesiveForce,
            SphericParticle *const element,
            Condition *const wall,
            bool &rSliding) override;

        double CalculateNormalForce(
            SphericParticle *const element1,
            SphericParticle *const element2,
            const double indentation,
            double LocalCoordSystem[3][3]) override;

        double CalculateNormalForce(
            SphericParticle *const element,
            Condition *const wall,
            const double indentation) override;

        double CalculateNormalForce(const double indentation) override;

        double CalculateCohesiveNormalForce(
            SphericParticle *const element1,
            SphericParticle *const element2,
            const double indentation) override;

        double CalculateCohesiveNormalForceWithFEM(
            SphericParticle *const element,
            Condition *const wall,
            const double indentation) override;

    private:
        /**
         * @brief Compact-support cubic-spline-like DPD kernel.
         *
         * r_surface is the distance measured from the DEM surface:
         * r_surface = center_distance - R_dem.
         *
         * For r_surface < 0, penetration is mapped to q = 0, giving the maximum
         * finite kernel value. This avoids a nonphysical kernel evaluation at
         * negative arguments while maintaining maximum conservative repulsion.
         */
        double KernelWeight(
            const double r_surface,
            const double h,
            const double rc) const;

        double DissipativeWeight(
            const double r_surface,
            const double h,
            const double rc) const;

        /**
         * @brief Returns the radius of the suspended DEM particle.
         *
         * The current simple implementation selects the larger particle radius.
         * This is correct when DPD radii are much smaller than suspended DEM radii.
         * Replace this by an explicit particle-type flag later if DPD particles can
         * ever be larger than the suspended DEM particle.
         */
        double GetDemRadius(
            SphericParticle *const element1,
            SphericParticle *const element2) const;

        double CalculateCenterDistance(
            SphericParticle *const element1,
            SphericParticle *const element2) const;

        friend class Serializer;

        void save(Serializer &rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(
                rSerializer,
                DEM_D_Linear_viscous_Coulomb);
        }

        void load(Serializer &rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(
                rSerializer,
                DEM_D_Linear_viscous_Coulomb);
        }
    };

} // namespace Kratos

#endif // KRATOS_DEM_DPD_DRAG_H_INCLUDED