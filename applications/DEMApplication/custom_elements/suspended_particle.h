#if !defined(KRATOS_SUSPENDED_PARTICLE_H_INCLUDED)
#define KRATOS_SUSPENDED_PARTICLE_H_INCLUDED

#include "custom_elements/spheric_particle.h"

namespace Kratos
{

class DPDParticle;

/**
 * @class SuspendedParticle
 * @ingroup DEMApplication
 * @brief DEM solid with DPD-only interaction against DPDParticle.
 *
 * Pair law:
 * - SuspendedParticle -- DPDParticle: DPD only, no base DEM contact.
 * - SuspendedParticle -- SuspendedParticle: ordinary DEM.
 * - SuspendedParticle -- SphericParticle: ordinary DEM.
 * - SuspendedParticle -- wall: ordinary DEM.
 */
class KRATOS_API(DEM_APPLICATION) SuspendedParticle
    : public SphericParticle
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SuspendedParticle);

    using SphericParticle::SphericParticle;

    SuspendedParticle() = default;

    SuspendedParticle(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : SphericParticle(NewId, pGeometry)
    {
    }

    SuspendedParticle(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : SphericParticle(NewId, pGeometry, pProperties)
    {
    }

    SuspendedParticle(
        IndexType NewId,
        NodesArrayType const& ThisNodes)
        : SphericParticle(NewId, ThisNodes)
    {
    }

    ~SuspendedParticle() override = default;

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<SuspendedParticle>(
            NewId,
            GetGeometry().Create(ThisNodes),
            pProperties);
    }

    void ComputeAdditionalForces(
        array_1d<double, 3>& rAdditionalForce,
        array_1d<double, 3>& rAdditionalMoment,
        const ProcessInfo& rProcessInfo,
        const array_1d<double, 3>& rGravity) override;

    std::string Info() const override
    {
        return "SuspendedParticle";
    }

protected:
    /**
     * @brief Returns false only for a DPDParticle neighbour.
     *
     * The base SphericParticle loop uses this method immediately before
     * its normal DEM contact-law/historical-force processing.
     */
    bool ShouldComputeDEMContactWith(
        const SphericParticle* pNeighbour) const override;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(
            rSerializer,
            SphericParticle);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(
            rSerializer,
            SphericParticle);
    }
};

} // namespace Kratos

#endif // KRATOS_SUSPENDED_PARTICLE_H_INCLUDED