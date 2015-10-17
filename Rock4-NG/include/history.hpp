#ifndef history__h
#define history__h
#include <utility>
#include <vector>
#include <iostream>
#include <string>
using namespace std;
namespace odes {
  ////////////////////////////////////////////////////////////////////////////
  /// Store the history of time steps: actually we store pairs <value,bool>
  /// where value is a time step used, and bool is tru iff it wa successfull.
  ///////////////////////////////////////////////////////////////////////////
  class history
  {
    static const int reserved=1000;
    typedef pair<double,bool> item;
    vector<item> V;
    //! help put an item in a human readable form.
    string sv(bool t)
    {
      if(t)
	return "accepted";
      else
	return "refused";
    }
  public:
    //! Constructor
    history()
    {
      V.reserve(reserved);
    }
    //! Destructor
    ~history(){}
    //! Enter a value
    //! \param h a time step.
    //! \param accepted (was the time step accepted?).
    void push(double h,bool accepted)
    {
      V.push_back(make_pair(h,accepted));
    }
    //! Reset, ie forget all.
    void reset()
    {
      V.clear(); V.reserve(reserved);
    }
    //! Return how many pairs are stored.
    int size() const{return V.size();}
    //! print.
    void print()
    {
      for(vector<item>::const_iterator I=V.begin();I!=V.end();I++)
	cout<<I->first<<" "<<sv(I->second)<<endl;
    }
    //! get the nth successfull time step, starting from the last
    // (ie: n===0 => last time step succesfull).
    double lastSuccess(int n)
    {
      vector<item>::const_reverse_iterator I=V.crbegin();
      int i=0; double ret;
      while(I!=V.crend() && i<=n)
	{
	  if(I->second)
	    {
	      ret=I->first; ++i;
	    }
	  ++I;
	}
      return ret;
    }
  };
};
#endif
