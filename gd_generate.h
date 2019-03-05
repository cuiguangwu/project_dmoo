#pragma once

#include "emo/IndividualMO.h"
#include "emo/PopulationMO.h"
#include "emo/Sel.h"
# include "global.h"
# include "global2.h"
# include "rand.h"

namespace az
{
	//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
	//!\brief namespace of dynamic evolutionary algoirhtm
namespace dea
{
class gd
{
public:
	void generater (unsigned int sizenew, az::mea::CPopulationMO& popnew, az::mea::CPopulationMO& pop);
};
}
}
}