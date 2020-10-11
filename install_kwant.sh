

if [ ! -f "/usr/local/bin/conda"  ]; then
echo "Installing packages for the activity"
export PYTHONPATH=""
MINICONDA_INSTALLER_SCRIPT=Miniconda3-latest-Linux-x86_64.sh
MINICONDA_PREFIX=/usr/local
wget https://repo.continuum.io/miniconda/$MINICONDA_INSTALLER_SCRIPT&>log
chmod +x $MINICONDA_INSTALLER_SCRIPT
./$MINICONDA_INSTALLER_SCRIPT -b -f -p $MINICONDA_PREFIX&>log
conda install -c conda-forge kwant &>log
conda install --channel defaults conda python=3.6 --yes&>log
conda update --channel defaults --all --yes&>log
fi
echo "All requested packages already installed."
