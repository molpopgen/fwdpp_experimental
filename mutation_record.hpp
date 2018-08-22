#ifndef FWDPP_ANCESTRY_MUTATION_RECORD_HPP
#define FWDPP_ANCESTRY_MUTATION_RECORD_HPP

#include <cstdint>

namespace fwdpp
{
	namespace ts
	{
		struct mutation_record
		{
			std::int32_t node;
			std::size_t key;
		};
	}
}

#endif
