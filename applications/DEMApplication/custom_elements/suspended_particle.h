#if !defined(KRATOS_SUSPENDED_PARTICLE_H_INCLUDED)
#define KRATOS_SUSPENDED_PARTICLE_H_INCLUDED

#include "custom_elements/spheric_particle.h"

namespace Kratos
{

/**
 * @class SuspendedParticle
 * @ingroup DEMApplication
 * @brief DEM solid particle with DPD interaction only against DPDParticle neighbours.
 *
 * Pair behavior:
 * - SuspendedParticle--DPDParticle:
 *     DPD range interaction, through a DEM_DPD_SPH_LIKE constitutive law.
 * - SuspendedParticle--SuspendedParticle:
 *     Ordinary SphericParticle DEM contact.
 * - SuspendedParticle--SphericParticle:
 *     Ordinary SphericParticle DEM contact.
 * - SuspendedParticle--rigid wall:
 *     Ordinary SphericParticle DEM wall contact.
 *
 * The DPD cross force is evaluated here, not by DPDParticle. Therefore
 * DPDParticle::ComputeBallToBallContactForceAndMoment must skip
 * SuspendedParticle neighbours.
 */
class KRATOS_API(DEM_APPLICATION) SuspendedParticle : public SphericParticle
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

    /**
     * @brief Hybrid particle--particle interaction loop.
     *
     * DPDParticle neighbours are handled with a DPD law over DPD_CUTOFF_RADIUS.
     * Every other neighbour is processed by the original SphericParticle
     * DEM contact implementation.
     */
    void ComputeBallToBallContactForceAndMoment(
        ParticleDataBuffer& rDataBuffer,
        const ProcessInfo& rProcessInfo,
        array_1d<double, 3>& rElasticForce,
        array_1d<double, 3>& rContactForce) override;

    /**
     * @brief Keep the normal DEM sphere--wall interaction unchanged.
     */
    void ComputeBallToRigidFaceContactForceAndMoment(
        ParticleDataBuffer& rDataBuffer,
        array_1d<double, 3>& rElasticForce,
        array_1d<double, 3>& rContactForce,
        array_1d<double, 3>& rRigidElementForce,
        const ProcessInfo& rProcessInfo) override;

    std::string Info() const override
    {
        return "SuspendedParticle";
    }

protected:
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle);
    }
};

} // namespace Kratos

#endif // KRATOS_SUSPENDED_PARTICLE_H_INCLUDED