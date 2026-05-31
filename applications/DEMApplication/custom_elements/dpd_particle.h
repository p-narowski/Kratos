#pragma once
#include "custom_elements/spheric_particle.h"
#include "DEM_application_variables.h"

namespace Kratos
{

    class KRATOS_API(DEM_APPLICATION) DPDParticle : public SphericParticle
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(DPDParticle);

        using SphericParticle::SphericParticle;

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const override
        {
            return Kratos::make_intrusive<DPDParticle>(NewId, GetGeometry().Create(ThisNodes), pProperties);
        }

        // Override the full ball-to-ball loop so DPD range forces fire for
        // ALL neighbours within r_cut, not only overlapping pairs.
        void ComputeBallToBallContactForceAndMoment(
            ParticleDataBuffer& data_buffer,
            const ProcessInfo& r_process_info,
            array_1d<double, 3>& rElasticForce,
            array_1d<double, 3>& rContactForce) override;

        std::string Info() const override { return "DPDParticle"; }

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
