#pragma once
#include "custom_elements/spheric_particle.h"
#include "DEM_application_variables.h"
#include <cmath>

namespace Kratos
{

class KRATOS_API(DEM_APPLICATION) SuspendedParticle : public SphericParticle
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SuspendedParticle);

    using SphericParticle::SphericParticle;
    using SphericParticle::ComputeAdditionalForces;

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<SuspendedParticle>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    std::string Info() const override { return "SuspendedParticle"; }

    void ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force,
                                 array_1d<double, 3>& additionally_applied_moment,
                                 const ProcessInfo& r_current_process_info,
                                 const array_1d<double,3>& gravity) override;

private:
    static double KernelWeight(const double r, const double h, const double rc)
    {
        if (r >= rc) return 0.0;
        const double q = r / h;
        if (q >= 0.0 && q < 1.0) {
            return 1.0 - 1.5 * q * q + 0.75 * q * q * q;
        } else if (q >= 1.0 && q < 2.0) {
            const double a = 2.0 - q;
            return 0.25 * a * a * a;
        }
        return 0.0;
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