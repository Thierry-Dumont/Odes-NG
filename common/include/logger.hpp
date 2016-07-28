#ifndef logger__h
#define logger__h
#include <list>
#include <iostream>
#include <string>
using namespace std;
namespace odes
{
  enum eventType{success=0, rejectedStep, NewtonFailed,
		 NewtonWillNotConverge,changedH,all};
  /////////////////////////////////////////////////////////////////////////
  /// For logging events (change of time steps, rejected steps, and so on).
  /// \brief For logging events.
  ////////////////////////////////////////////////////////////////////////
  class logger
  {
  private:
    struct event
    {
      double time, value;
      int step;
      eventType TheEvent;
      event(double _time,int _step,eventType _TheEvent,double _value=0.0)
      {
    	time=_time; step=_step; TheEvent=_TheEvent;value=_value;
      }
    };
    list<event> L;
    
  public:
    //! constructor
    logger(){}
    //! destructor
    ~logger(){}
    //! reset.
    void clear() 
    {
      L.clear();
    }
    //! put an event
    //!\param _time actual time
    //!\param _step number
    //!\param _TheEvent event type (see struct event).
    //!\param _value some value like new time step.
    void put(double _time,int _step,eventType _TheEvent,double _value)
    {
      L.push_back(event(_time,_step,_TheEvent, _value));
    }
    //! convert to a string.
    //!\param E filter on a given  type of event; if E==all, use all events. 
    string tostring(eventType E=all) const
    {
      string names[6]={"success","rejectedStep","NewtonFailed",
			      "NewtonWillNotConverge","changedH","all"};
      string ret="";
      int count=0;
      for(list<event>::const_iterator I=L.begin();I!=L.end();I++)
	{ 
	  event Ev=*I;
	  if(Ev.TheEvent==E|| E==all)
	    {
	      string tmp=to_string(++count)+
	       	" Time: "+to_string(Ev.time)+" step: "+to_string(Ev.step)
	       	+" Event: "+names[Ev.TheEvent];
	      if(Ev.TheEvent==changedH)
	      	tmp+=" new step: "+to_string(Ev.value);
	      tmp+="\n";
	      ret+=tmp;
	    }
	}
      return ret;
    }
    //! print
    //! \param E eventype, for filtering results.
    void print(eventType E=all)
    {
 
      string s=tostring(all);
      if(s.size()>0)
	cout<<s;
      else
	cout<<"No event."<<endl;
    }
  };
  ///! overload <<
  std::ostream& operator<<(std::ostream& os, const logger& obj)
  {
    // write obj to stream
    os<<obj.tostring();
    return os;
  }
}
#endif
