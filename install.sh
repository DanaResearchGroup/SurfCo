#!/bin/bash

# SurfCo - Complete Installation Script
# Protein Corona Prediction Framework
# Installs Python dependencies and UnitedAtom component

echo "=================================================================="
echo "   SurfCo Installation Script"
echo "   Protein Corona Prediction Framework"
echo "   Complete setup for all components"
echo "=================================================================="
echo ""

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to detect OS
detect_os() {
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        echo "linux"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        echo "macos"
    else
        echo "unsupported"
    fi
}

# Function to install system dependencies
install_system_deps() {
    local os="$1"
    echo "Installing system dependencies..."
    
    if [ "$os" == "linux" ]; then
        # Detect Linux distribution
        if command_exists apt-get; then
            echo "  Detected Debian/Ubuntu system"
            echo "  Installing: build-essential libboost-all-dev"
            sudo apt-get update
            sudo apt-get install -y build-essential libboost-all-dev
        elif command_exists yum; then
            echo "  Detected RHEL/CentOS system"
            echo "  Installing: gcc-c++ make boost-devel"
            sudo yum install -y gcc-c++ make boost-devel
        elif command_exists dnf; then
            echo "  Detected Fedora system"
            echo "  Installing: gcc-c++ make boost-devel"
            sudo dnf install -y gcc-c++ make boost-devel
        else
            echo "  Unknown Linux distribution. Please install:"
            echo "    - C++ compiler (g++)"
            echo "    - Make build system"
            echo "    - Boost C++ libraries"
            return 1
        fi
    elif [ "$os" == "macos" ]; then
        if command_exists brew; then
            echo "  Installing Boost via Homebrew..."
            brew install boost
        else
            echo "  Homebrew not found. Installing Homebrew first..."
            /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
            brew install boost
        fi
    fi
    
    return 0
}

# Start installation
OS=$(detect_os)
echo "Detected OS: $OS"
echo ""

# Check if we're in the SurfCo directory
if [ ! -f "run_pipeline.py" ]; then
    echo "❌ Error: Please run this script from the SurfCo root directory"
    echo "   (The directory containing run_pipeline.py)"
    exit 1
fi

echo "Step 1: Checking Python installation..."
echo "----------------------------------------"
if command_exists python3; then
    PYTHON_VERSION=$(python3 --version 2>&1 | grep -Po '(?<=Python )\d+\.\d+')
    echo "✅ Python $PYTHON_VERSION found"
else
    echo "❌ Python 3 not found. Please install Python 3.7 or higher"
    exit 1
fi

echo ""
echo "Step 2: Setting up Python environment..."
echo "----------------------------------------"
read -p "Create virtual environment? (recommended) [y/n]: " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    if ! command_exists virtualenv && ! python3 -m venv --help &> /dev/null; then
        echo "Installing virtualenv..."
        pip3 install virtualenv
    fi
    
    echo "Creating virtual environment..."
    python3 -m venv surfco_env
    
    echo "Activating virtual environment..."
    source surfco_env/bin/activate
    
    echo "✅ Virtual environment created and activated"
    VENV_CREATED=true
else
    VENV_CREATED=false
fi

echo ""
echo "Step 3: Installing Python dependencies..."
echo "----------------------------------------"
pip install --upgrade pip
if pip install -r requirements.txt; then
    echo "✅ Python dependencies installed successfully"
else
    echo "❌ Failed to install Python dependencies"
    exit 1
fi

echo ""
echo "Step 4: Checking system dependencies for UnitedAtom..."
echo "-----------------------------------------------------"

MISSING_DEPS=""

# Check for C++ compiler
if ! command_exists g++; then
    MISSING_DEPS="$MISSING_DEPS g++"
fi

# Check for make
if ! command_exists make; then
    MISSING_DEPS="$MISSING_DEPS make"
fi

# Check for boost libraries (Linux)
if [ "$OS" == "linux" ]; then
    if ! ldconfig -p | grep -q libboost_filesystem 2>/dev/null; then
        MISSING_DEPS="$MISSING_DEPS libboost-dev"
    fi
elif [ "$OS" == "macos" ]; then
    if ! command_exists brew; then
        echo "⚠️  Homebrew not found. Will attempt to install later."
        MISSING_DEPS="$MISSING_DEPS boost"
    elif ! brew list boost &>/dev/null; then
        MISSING_DEPS="$MISSING_DEPS boost"
    fi
fi

if [ -n "$MISSING_DEPS" ]; then
    echo "⚠️  Missing system dependencies:$MISSING_DEPS"
    echo ""
    read -p "Attempt to install system dependencies? [y/n]: " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        install_system_deps "$OS"
        if [ $? -ne 0 ]; then
            echo "❌ Failed to install system dependencies"
            echo "   Please install manually and re-run this script"
            exit 1
        fi
        echo "✅ System dependencies installed"
    else
        echo "⚠️  Continuing without installing system dependencies"
        echo "   UnitedAtom compilation may fail"
    fi
else
    echo "✅ All system dependencies found"
fi

echo ""
echo "Step 5: Installing UnitedAtom..."
echo "--------------------------------"

# Check if UA directory exists
if [ ! -d "ua" ]; then
    echo "❌ UnitedAtom source directory 'ua' not found"
    echo "   Please ensure the 'ua' directory with UnitedAtom source is present"
    exit 1
fi

# Change to UA directory
cd ua

# Check if UnitedAtom binary already exists
if [ -f "build/UnitedAtom" ] || [ -f "UnitedAtom" ]; then
    echo "UnitedAtom binary found. Skipping compilation."
    UA_EXISTS=true
else
    echo "Compiling UnitedAtom..."
    
    # Clean any previous builds
    make clean 2>/dev/null || true
    rm -f UnitedAtom *.o 2>/dev/null || true
    
    # Create build directory if it doesn't exist
    mkdir -p build
    
    # Try to compile
    echo "  Building UnitedAtom (this may take a few minutes)..."
    if make > build.log 2>&1; then
        if [ -f "build/UnitedAtom" ]; then
            echo "✅ UnitedAtom compiled successfully!"
            UA_EXISTS=true
        elif [ -f "UnitedAtom" ]; then
            # Move to build directory for consistency
            mv UnitedAtom build/
            echo "✅ UnitedAtom compiled successfully!"
            UA_EXISTS=true
        else
            echo "❌ UnitedAtom compilation completed but binary not found"
            echo "   Check build.log for details"
            UA_EXISTS=false
        fi
    else
        echo "❌ UnitedAtom compilation failed"
        echo "   Error details:"
        tail -10 build.log
        echo ""
        echo "   This may be due to:"
        echo "   1. Missing C++ compiler or development tools"
        echo "   2. Missing Boost libraries"
        echo "   3. Incompatible source code"
        echo ""
        echo "   You can:"
        echo "   1. Install missing dependencies and re-run this script"
        echo "   2. Continue installation (SurfCo will work without energy calculations)"
        UA_EXISTS=false
    fi
fi

# Return to root directory
cd ..

# Create symbolic link if UA was built successfully
if [ "$UA_EXISTS" = true ]; then
    if [ -f "ua/build/UnitedAtom" ]; then
        ln -sf ua/build/UnitedAtom UnitedAtom
        echo "✅ Created UnitedAtom symbolic link"
    elif [ -f "ua/UnitedAtom" ]; then
        ln -sf ua/UnitedAtom UnitedAtom
        echo "✅ Created UnitedAtom symbolic link"
    fi
fi

echo ""
echo "Step 6: Setting up directory structure..."
echo "----------------------------------------"
mkdir -p data input projects

# Check for essential data files
if [ ! -f "data/MaterialSet.csv" ]; then
    echo "⚠️  Warning: data/MaterialSet.csv not found"
    echo "   Please ensure you have the necessary PMF and Hamaker data files"
    echo "   in the data/ directory for UnitedAtom to function properly"
fi

echo "✅ Directory structure created"

echo ""
echo "Step 7: Testing installation..."
echo "------------------------------"

# Test Python imports
echo "Testing Python modules..."
python3 -c "
import sys
try:
    import numpy, pandas, scipy, matplotlib, yaml
    print('✅ All Python modules imported successfully')
except ImportError as e:
    print(f'❌ Python module import failed: {e}')
    sys.exit(1)
"

if [ $? -ne 0 ]; then
    echo "❌ Python module test failed"
    exit 1
fi

# Test UnitedAtom if available
if [ "$UA_EXISTS" = true ] && [ -f "UnitedAtom" ]; then
    echo "Testing UnitedAtom..."
    if ./UnitedAtom --help >/dev/null 2>&1 || ./UnitedAtom >/dev/null 2>&1; then
        echo "✅ UnitedAtom executable works"
    else
        echo "⚠️  UnitedAtom executable found but may have issues"
    fi
else
    echo "⚠️  UnitedAtom not available - energy calculations will be disabled"
fi

# Final summary
echo ""
echo "=================================================================="
echo "   Installation Complete!"
echo "=================================================================="
echo ""

if [ "$VENV_CREATED" = true ]; then
    echo "📋 To run SurfCo:"
    echo "   1. Activate virtual environment: source surfco_env/bin/activate"
    echo "   2. Place your PDB files and config.yaml in: input/[project_name]/"
    echo "   3. Run: python run_pipeline.py [project_name]"
    echo ""
    echo "💡 Remember to activate the virtual environment before each use:"
    echo "   source surfco_env/bin/activate"
else
    echo "📋 To run SurfCo:"
    echo "   1. Place your PDB files and config.yaml in: input/[project_name]/"
    echo "   2. Run: python run_pipeline.py [project_name]"
fi

echo ""
echo "📚 For help and documentation:"
echo "   python run_pipeline.py --help"
echo "   python run_pipeline.py --show-c    # Citation requirements"
echo ""

if [ "$UA_EXISTS" = true ]; then
    echo "✅ Full installation successful - all components available"
else
    echo "⚠️  Partial installation - UnitedAtom not available"
    echo "   SurfCo will work but energy calculations will be disabled"
    echo "   You can re-run this script after installing system dependencies"
fi

echo ""
echo "🔗 Documentation: https://github.com/DanaResearchGroup/SurfCo"
echo "=================================================================="