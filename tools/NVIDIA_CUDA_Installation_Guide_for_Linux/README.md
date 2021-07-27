# NVIDIA CUDA Installation Guide for Linux Ubuntu

**Dated version: CUDA Toolkit 11.4**

This summary is specifically prepared for Ubuntu workstations or laptops equipped with NVIDIA CUDA-capable GPU(s).

The procedures for all supported Linux distros are on [https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html).

## Pre-installation Actions

### Install build essentials

	$ sudo apt-get install build-essential ssh libopenmpi-dev

### Keep the system up-to-date

	$ sudo apt update
	$ sudo apt upgrade
	$ sudo apt dist-upgrade

### Verify You Have a CUDA-Capable GPU

	$ lspci | grep -i nvidia

### Download the NVIDIA CUDA Toolkit

	$ cd ~/Downloads
	$ wget https://developer.download.nvidia.com/compute/cuda/11.4.0/local_installers/cuda_11.4.0_470.42.01_linux.run

### Handle Conflicting Installation Methods

Before installing CUDA, any previously installations that could conflict should be uninstalled. This will not affect systems which have not had CUDA installed previously, or systems where the installation method has been preserved (RPM/Deb vs. Runfile).

Use the following command to uninstall a Toolkit runfile installation:

	$ sudo /usr/local/cuda-X.Y/bin/cuda-uninstaller

Use the following command to uninstall a Driver runfile installation:
 
	$ sudo /usr/bin/nvidia-uninstall
 
Use the following commands to uninstall a RPM/Deb installation:
 
	$ sudo apt-get --purge remove <package_name>
 
## Runfile Installation

The Runfile installation installs the NVIDIA Driver, CUDA Toolkit, and CUDA Samples via an interactive ncurses-based interface.

### Disable the Nouveau drivers

To install the Display Driver, the Nouveau drivers must first be disabled.

The Nouveau drivers are loaded if the following command prints anything:

	$ lsmod | grep nouveau

Create a file at `/etc/modprobe.d/blacklist-nouveau.conf`

	sudo vi /etc/modprobe.d/blacklist-nouveau.conf

with the following contents:

	blacklist nouveau
	options nouveau modeset=0
	
Regenerate the kernel initramfs:

	$ sudo update-initramfs -u
	
### Reboot into text mode (runlevel 3).

This can usually be accomplished by adding the number "3" to the end of the system's kernel boot parameters.

Since the NVIDIA drivers are not yet installed, the text terminals may not display correctly. Temporarily adding "nomodeset" to the system's kernel boot parameters may fix this issue.

The reboot is required to completely unload the Nouveau drivers and prevent the graphical interface from loading. The CUDA driver cannot be installed while the Nouveau drivers are loaded or while the graphical interface is active.

**Important**: Verify that the Nouveau drivers are not loaded. It should show nothing with the following command.

	$ lsmod | grep nouveau
	
Run the installer and follow the on-screen prompts:

	$ cd ~/
	$ sudo sh ./cuda_<version>_linux.run # example: cuda_11.4.0_470.42.01_linux.run

The installer will prompt for the following:

* EULA Acceptance
* CUDA Driver installation
* CUDA Toolkit installation, location, and `/usr/local/cuda` symbolic link
* CUDA Samples installation and location

The default installation locations for the toolkit and samples are:

* **CUDA Toolkit**: `/usr/local/cuda-11.4`
* **CUDA Samples**: `$(HOME)/NVIDIA_CUDA-11.4_Samples`

**Note**: Installing Mesa may overwrite the /usr/lib/libGL.so that was previously installed by the NVIDIA driver, so a reinstallation of the NVIDIA driver might be required after installing these libraries.

Reboot the system to reload the graphical interface.

## Post-installation Actions

### Environment Setup

To add this path to the PATH variable (e.g., in `~/.bashrc` file):

	$ export PATH=/usr/local/cuda-11.4/bin${PATH:+:${PATH}}
	$ export LD_LIBRARY_PATH=/usr/local/cuda-11.4/lib64\
                         ${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

### Recommended Actions

The version of the CUDA Toolkit can be checked by running

	$ nvcc -V

in a terminal window. The `nvcc` command runs the compiler driver that compiles CUDA programs. It calls the gcc compiler for C code and the NVIDIA PTX compiler for the CUDA code.

The NVIDIA CUDA Toolkit includes sample programs in source form. You should compile them by changing to

	$ cd ~/NVIDIA_CUDA-11.4_Samples

and typing

	$ make
	
The resulting binaries will be placed under `~/NVIDIA_CUDA-11.4_Samples/bin`.

After compilation, find and run `deviceQuery` under `~/NVIDIA_CUDA-11.4_Samples`. If the CUDA software is installed and configured correctly, the output for `deviceQuery` should look similar to that shown in

[https://docs.nvidia.com/cuda/cuda-installation-guide-linux/graphics/valid-results-from-sample-cuda-devicequery-program.png](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/graphics/valid-results-from-sample-cuda-devicequery-program.png)