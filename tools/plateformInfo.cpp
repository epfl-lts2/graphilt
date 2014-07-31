/****************************************************************************
**
** Copyright (C) 2013 EPFL-LTS2
** Contact: Kirell Benzi (first.last@epfl.ch)
**
** This file is part of Graphilt.
**
**
** GNU General Public License Usage
** This file may be used under the terms of the GNU
** General Public License version 3.0 as published by the Free Software
** Foundation and appearing in the file LICENSE.GPLv3 included in the
** packaging of this file.  Please review the following information to
** ensure the GNU General Public License version 3.0 requirements
** will be met: http://www.gnu.org/licenses/
**
****************************************************************************/

#include <iostream>

#include <viennacl/ocl/device.hpp>
#include <viennacl/ocl/platform.hpp>

int main()
{
   typedef std::vector< viennacl::ocl::platform > platforms_type;
   platforms_type platforms = viennacl::ocl::get_platforms();
   
   bool is_first_element = true;
   for( platforms_type::iterator platform_iter  = platforms.begin();
                                 platform_iter != platforms.end();
                               ++platform_iter )
   {
    typedef std::vector<viennacl::ocl::device> devices_type;
    devices_type devices = platform_iter->devices(CL_DEVICE_TYPE_ALL);
    
    //
    // print some platform info
    //
    std::cout << "# =========================================" << std::endl;
    std::cout << "#         Platform Information             " << std::endl;
    std::cout << "# =========================================" << std::endl;
    
    std::cout << "#" << std::endl;
    std::cout << "# Vendor and version: " << platform_iter->info() << std::endl;
    
    if( is_first_element ) {
      std::cout << "# ViennaCL uses this OpenCL platform by default." << std::endl;
      is_first_element = false;
    }

    std::cout << "# " << std::endl;
    std::cout << "# Available Devices: " << std::endl;
    for( devices_type::iterator iter = devices.begin(); iter != devices.end(); iter++ ) {
        std::cout << std::endl;

        std::cout << "  -----------------------------------------" << std::endl;
        std::cout << "  No.:              " << std::distance(devices.begin(), iter) << std::endl;
        std::cout << "  Name:             " << iter->name() << std::endl;
        std::cout << "  Compute Units:    " << iter->max_compute_units() << std::endl;
        std::cout << "  Workgroup Size:   " << iter->max_workgroup_size() << std::endl;
        std::cout << "  Global Memory:    " << iter->global_memory()/(1024*1024) << " MB" << std::endl;
        std::cout << "  Local Memory:     " << iter->local_memory()/1024 << " KB" << std::endl;
        std::cout << "  Max-alloc Memory: " << iter->max_allocable_memory()/(1024*1024) << " MB" << std::endl;
        std::cout << "  Double Support:   " << iter->double_support() << std::endl;
        std::cout << "  Driver Version:   " << iter->driver_version() << std::endl;
        std::cout << "  -----------------------------------------" << std::endl;
    }
    std::cout << std::endl;
    std::cout << "###########################################" << std::endl;
    std::cout << std::endl;
   }
   
   return EXIT_SUCCESS;
}



