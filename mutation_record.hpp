#ifndef FWDPP_ANCESTRY_MUTATION_RECORD_HPP__
#define FWDPP_ANCESTRY_MUTATION_RECORD_HPP__

#include <cstdint>

namespace fwdpp
{
	namespace ancestry
	{
		struct mutation_record
		{
			std::int32_t node;
			std::size_t key;
		};
	}
}

#endif
