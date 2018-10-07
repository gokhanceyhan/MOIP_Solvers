//
//  subspace.h
//  FirstProject
//
//  Created by Gokhan Ceyhan on 3/25/15.
//  Copyright (c) 2015 Gokhan Ceyhan. All rights reserved.
//

#ifndef subspace_h
#define subspace_h

#include <vector>


class subspace {
    
	private:

	std::vector<float> bound;
    std::vector<float> point;
    int sub_status; // -1: initial, 0: infeasible, 1: feasible
    int mod_status; // -1: initial, 0: no model solved, 1: a model is solved to obtain sub_status info

    public:
    
    inline void set_bound(std::vector<float>& b);
    
    inline std::vector<float> get_bound ();
    
    inline void set_point(std::vector<float>& p);
    
    inline std::vector<float> get_point ();
    
    inline void set_sub_status(int s_status);
    
    inline int get_sub_status();
    
    inline void set_mod_status(int m_status);
    
    inline int get_mod_status();

  
};

inline void subspace::set_bound(std::vector<float>& b)
{
    bound=b;
}

inline std::vector<float> subspace::get_bound ()
{
    return bound;
}

inline void subspace::set_point(std::vector<float>& p)
{
    point=p;
}

inline std::vector<float> subspace::get_point ()
{
    return point;
}

inline void subspace::set_sub_status(int s_status)
{
    sub_status =s_status;
}

inline int subspace::get_sub_status()
{
    return sub_status;
}

inline void subspace::set_mod_status(int m_status)
{
    mod_status =m_status;
}

inline int subspace::get_mod_status()
{
    return mod_status;
}



#endif
