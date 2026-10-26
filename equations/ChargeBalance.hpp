#pragma once

#include "equations/Equation.hpp"

class ChargeBalance : public Equation
{
protected:
  ParFiniteElementSpace & l2_fespace;
  const Region particle_region;
  ParMixedBilinearForm * C = nullptr;
  HypreParMatrix Cmat;

public:
  ChargeBalance(ParFiniteElementSpace & h1, ParFiniteElementSpace & l2, const Region & region)
    : Equation(h1), l2_fespace(l2), particle_region(region)
  {
    
   if(particle_region == NE)
   {
    Array<int> ess_bdr(fespace.GetParMesh()->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[3 - 1] = 1; 
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

 }
  

  virtual void Update();
  const HypreParMatrix & GetC() const { return Cmat; };

  ~ChargeBalance() override;
};
