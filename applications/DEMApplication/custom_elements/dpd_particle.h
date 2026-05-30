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

        Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes,
                                PropertiesType::Pointer pProperties) const override
        {
            return Kratos::make_intrusive<DPDParticle>(NewId, GetGeometry().Create(ThisNodes), pProperties);
        }

        // Returns rc for wall/neighbour search — makes law radius-agnostic
        double GetInteractionRadius(const int radius_index = 0) override
        {
            if (GetProperties().Has(DPD_CUTOFF_RADIUS))
                return GetProperties()[DPD_CUTOFF_RADIUS];
            return SphericParticle::GetInteractionRadius(radius_index);
        }

        std::string Info() const override { return "DPDParticle"; }

    protected:
        friend class Serializer;
        void save(Serializer &rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle);
        }
        void load(Serializer &rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle);
        }
    };

} // namespace Kratos