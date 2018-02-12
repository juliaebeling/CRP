#ifndef CRPROPA_HADRONICINTERACTION_H
#define CRPROPA_HADRONICINTERACTION_H

#include "crpropa/Module.h"

#include <vector>

namespace crpropa {
    
    
    class HadronicInteraction: public Module {
        
        
        public:
        HadronicInteraction( );
        void process(Candidate *candidate) const;
        double distribution_my1(double energy, double x) const;
        double distribution_e(double energy, double x) const;
        double distribution_g(double energy, double x) const;
        //~ double counter(Candidate *candidate, double density) const;
        //~ double counterPion(Candidate *candidate, double density) const;
    };
    
} // namespace crpropa

#endif // CRPROPA_HADRONICINTERACTION_H
