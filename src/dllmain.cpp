// dllmain.cpp : Defines the entry point for the DLL application.
#include "stdafx.h"
// This file defines the required functions for registering the plugin classes
// with the FEBio framework. 
#include "FECore/FECoreKernel.h"
#include "Holzapfel-Ogden.h"

//-----------------------------------------------------------------------------
// This required function returns the version of the FEBio SDK that is being
// used by the plugin. This version number will be checked by FEBio to make
// sure that the FEBio executable is compatible with this SDK. Usually this
// function just returns the predefined macro FE_SDK_VERSION.
FECORE_EXPORT unsigned int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

//-----------------------------------------------------------------------------
// If FEBio is successful in loading the plugin, this is the first function
// that will be called and can be used to initialize any resources that are
// needed by the plugin. FEBio passes a reference to the FECoreKernel as a parameter,
// which can be used to access some of the FEBio resources such as the log file.
// This function is optional.
FECORE_EXPORT void PluginInitialize(FECoreKernel& febio)
{
	FECoreKernel::SetInstance(&febio);
    
    // Register the class. This macro requires the name of the class,
	// the super class ID, and a name for the class that will be used
	// in the FEBio input file. 
	REGISTER_FECORE_CLASS(HolzapfelMyocardiumPI, "Holzapfel_Ogden");
}

//-----------------------------------------------------------------------------
// This required function should return the number of classes that this plugin defines.
// This number determines the number of times that FEBio calls the RegisterPlugin
// later, which is used to register the actual factory classes.
//FECORE_EXPORT int PluginNumClasses()
//{
	//return 1;
//}

//-----------------------------------------------------------------------------
// This required function is called by FEBio and will be called the same number of times
// as return by PluginNumClasses. The parameter i is used as a counter which contains
// how often this function has been called thus far.
//FECORE_EXPORT FECoreFactory* PluginGetFactory(int i)
//{
	//if (i == 0) return &mat_factory;
	//return 0;
//}

//-----------------------------------------------------------------------------
// This function is called when FEBio exits and gives the plugin a chance to close
// any resource that have been allocated. 
// This is an optional function.
FECORE_EXPORT void PluginCleanup()
{

}
