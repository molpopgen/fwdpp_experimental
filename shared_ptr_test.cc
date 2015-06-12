#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

int main(int argc, char ** argv)
{
  using ptr_t = shared_ptr<int>;
  vector< ptr_t > mvector;

  for(int i = 0 ; i < 100 ; ++i )
    {
      mvector.push_back( ptr_t(new int(i)) );
    }

  //cool, now copy into "gametes"
  vector<ptr_t> mvector2;
  for( const auto & ptr : mvector ) mvector2.push_back(ptr);

  for(unsigned i=0;i<mvector.size();++i)
    {
      cout << mvector[i].use_count() << ' ' << mvector2[i].use_count() << ' '
	   << *mvector[i] << ' ' << *mvector2[i] << '\n';
    }

  //What happens upon sorting mvector by ascending order?
  sort(mvector.begin(),mvector.end(),[](const ptr_t & __p1,const ptr_t & __p2) { return *__p1 > *__p2; });
  for(unsigned i=0;i<mvector.size();++i)
    {
      cout << mvector[i].use_count() << ' ' << mvector2[i].use_count() << ' '
	   << *mvector[i] << ' ' << *mvector2[i] << '\n';
    }

  //Now, for the kicker.  Let's delete 1/2 the stuff in mvector
  mvector.erase( remove_if(mvector.begin(),mvector.end(),[](const ptr_t & __p) { return *__p <= 50; }), mvector.end() );

  for( unsigned i = 0 ; i< mvector2.size() ; ++i )
    {
      cout << mvector2[i].use_count() << ' ' << *mvector2[i] << '\n';
    }
}
