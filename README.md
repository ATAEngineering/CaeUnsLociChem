# Loci/CHEM Plugin for Pointwise
This plugin was developed by [ATA Engineering](http://www.ata-e.com) to 
export grids from Pointwise directly to Loci/CHEM's `vog` format. Using CAE
export, portions of the `vars` file may also be generated. It should be noted 
that this plugin does not optimize the grid layout as is done by the Loci
grid converstion utilities (e.g. `ccm2vog`, `ugrid2vog`, etc). Testing has
indicated that this lack of optimization results in roughly a 10% slow down
compared to an optimized grid.

# Installation
Building the plugin requires the [Poinwise Plugin SDK](https://www.pointwise.com/plugins/html/index.html).
The plugin also requres [HDF5](https://www.hdfgroup.org/solutions/hdf5/), it was tested 
with version 1.12.2, but other versions may work as well.

The SDK must first be used to create a plugin project for Loci/CHEM. 
See [here](https://www.pointwise.com/plugins/html/dd/d8c/create_cae_plugin_project.html) for more details,
but in general plugin creation from the `PluginSDK` directory should look like below:

```bash
./mkplugin -caeu -cpp LociChem
```

Next the source code from this repository must be integrated with the plugin code generated 
by the SDK in the `src/plugins/CaeUnsLociChem` directory. See 
[here](https://github.com/pointwise/How-To-Integrate-Plugin-Code) for more information. The make
files generated by the SDK will likely have to be edited to include the path to HDF5's `include`
and `lib` directories. The plugin should be linked against the following HDF5 libraries `libhdf5`,
`libz`, `libszaec`, and `libaec`.

Finally, the plugin must be built and the resulting library placed in the `PWI_PLUGINS_SEARCH_PATH`.
See [here](https://www.pointwise.com/plugins/html/da/dde/build_cae_plugin.html) for more information.

# Updating to a new SDK
To update the plugin to a new SDK, which is occasionally required to support the latest versions
of Pointwise, follow these instructions for Windows.

* Extract the new SDK
* Update `PluginSDK/src/plugins/site.h` with group name and ID.
* Copy `CaeUnsLociChem` folder into `PLuginSDK/src/plugins/`
* Open Visual Studio `sln` file in `PluginSDK` directory
* Right click on `Plugins` in source tree, select `Add`, select `Existing Project`
* Pick `vcxproj` file inside `PLuginSDK/src/plugins/CaeUnsLociChem`
* To upgrade version of Visual Studio, right click on imported plugin name, select `Retarget` 