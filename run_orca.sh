#!/bin/bash
# =============================================================================
#                    ORCA 6.1 Workflow Script
# =============================================================================



# =============================================================================
#                         USER CONFIGURATION
# =============================================================================
# Modify these variables according to your environment

# ORCA installation path (REQUIRED)
ORCA_PATH="/opt/orca/orca_6_1_1"

# Default number of CPU cores per job (can be changed interactively)
DEFAULT_NPROCS=6

# Memory per core (MB) - Total memory = NPROCS * MAXCORE
MAXCORE=4000

# Default solvent for CPCM (DMF, THF, Water, Acetonitrile, or empty for gas phase)
DEFAULT_SOLVENT="DMF"

# MO Cube file settings
CUBE_GRID_X=80
CUBE_GRID_Y=80
CUBE_GRID_Z=80
DEFAULT_MO_RANGE=3    # HOMO-N to LUMO+N

# SLURM settings
SLURM_PARTITION=""           # Leave empty for default partition
SLURM_TIME="120:00:00"        # Max walltime
SLURM_MEM="100G"              # Total memory request (auto-calculated based on NPROCS)

# =============================================================================
#            DO NOT MODIFY BELOW UNLESS YOU KNOW WHAT YOU'RE DOING
# =============================================================================

# Colors - Text
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
NC='\033[0m'

# Colors - Background
BG_GREEN='\033[42m'
BG_RED='\033[41m'
BG_YELLOW='\033[43m'
BG_MAGENTA='\033[45m'
BG_BLUE='\033[44m'

# =============================================================================
# Function: Check SLURM cluster resources
# =============================================================================
check_slurm_resources() {
    echo -e "${BG_BLUE}${WHITE}====================================${NC}"
    echo -e "${BG_BLUE}${WHITE}     SLURM Cluster Resources        ${NC}"
    echo -e "${BG_BLUE}${WHITE}====================================${NC}"

    # Check if sinfo command is available
    if ! command -v sinfo &> /dev/null; then
        echo -e "${YELLOW}[WARN]${NC} sinfo command not found. Using default cores: ${DEFAULT_NPROCS}"
        NPROCS=${DEFAULT_NPROCS}
        return
    fi

    # Get cluster resource information
    echo -e "${CYAN}Partition Status:${NC}"
    echo "------------------------------------------------------------------------"

    # Show partition summary (PARTITION, AVAIL, NODES, STATE, CPUS)
    sinfo -h -o "%12P %5a %5D %10T %5c" 2>/dev/null | head -10

    echo "------------------------------------------------------------------------"

    # Calculate total and available cores
    TOTAL_CORES=$(sinfo -h -o "%C" 2>/dev/null | head -1)
    if [ -n "$TOTAL_CORES" ]; then
        # Format: ALLOCATED/IDLE/OTHER/TOTAL
        ALLOC=$(echo "$TOTAL_CORES" | cut -d'/' -f1)
        IDLE=$(echo "$TOTAL_CORES" | cut -d'/' -f2)
        OTHER=$(echo "$TOTAL_CORES" | cut -d'/' -f3)
        TOTAL=$(echo "$TOTAL_CORES" | cut -d'/' -f4)

        echo -e "Total Cores:     ${CYAN}${TOTAL}${NC}"
        echo -e "Allocated:       ${RED}${ALLOC}${NC}"
        echo -e "Idle (Free):     ${GREEN}${IDLE}${NC}"
        echo -e "Other:           ${YELLOW}${OTHER}${NC}"
    else
        IDLE="N/A"
        echo -e "${YELLOW}[WARN]${NC} Could not retrieve cluster information"
    fi

    # Show user's current jobs
    echo ""
    echo -e "${CYAN}Your Running Jobs:${NC}"
    USER_JOBS=$(squeue -u $USER -h -o "%i %j %t %C %M" 2>/dev/null | head -5)
    if [ -n "$USER_JOBS" ]; then
        echo "JOBID      NAME          ST  CPUS  TIME"
        echo "$USER_JOBS"
    else
        echo "(No running jobs)"
    fi

    echo "------------------------------------------------------------------------"

    # Prompt for number of cores
    echo ""
    if [ "$IDLE" != "N/A" ] && [ "$IDLE" -gt 0 ] 2>/dev/null; then
        echo -e "Available cores: ${GREEN}${IDLE}${NC}"
    fi
    echo -e "Default cores: ${CYAN}${DEFAULT_NPROCS}${NC}"
    echo ""

    while true; do
        read -p "Number of cores to use [default: ${DEFAULT_NPROCS}]: " NPROCS_INPUT
        NPROCS=${NPROCS_INPUT:-$DEFAULT_NPROCS}

        # Validate input is a positive integer
        if [[ "$NPROCS" =~ ^[0-9]+$ ]] && [ "$NPROCS" -gt 0 ]; then
            # Warn if requesting more than available
            if [ "$IDLE" != "N/A" ] && [ "$IDLE" -gt 0 ] 2>/dev/null; then
                if [ "$NPROCS" -gt "$IDLE" ]; then
                    echo -e "${BG_YELLOW}${WHITE} WARNING: Requesting ${NPROCS} cores, but only ${IDLE} idle. Job may be queued. ${NC}"
                    read -p "Continue anyway? [y/n]: " CONTINUE
                    if [[ ! "$CONTINUE" =~ ^[Yy]$ ]]; then
                        continue
                    fi
                fi
            fi
            break
        else
            echo -e "${BG_RED}${WHITE} Invalid input. Please enter a positive number. ${NC}"
        fi
    done

    echo -e "Selected cores: ${GREEN}${NPROCS}${NC}"

    # Update SLURM memory based on cores
    SLURM_MEM="$((NPROCS * MAXCORE / 1000 + 2))G"
    echo -e "Total memory: ${CYAN}${SLURM_MEM}${NC} (${MAXCORE}MB/core + 2GB buffer)"
}

# =============================================================================
# Function: Calculate total electrons from xyz file (pure bash)
# =============================================================================
get_total_electrons() {
    local xyz_file=$1

    declare -A ATOMIC_NUMBERS=(
        [H]=1 [He]=2 [Li]=3 [Be]=4 [B]=5 [C]=6 [N]=7 [O]=8 [F]=9 [Ne]=10
        [Na]=11 [Mg]=12 [Al]=13 [Si]=14 [P]=15 [S]=16 [Cl]=17 [Ar]=18
        [K]=19 [Ca]=20 [Sc]=21 [Ti]=22 [V]=23 [Cr]=24 [Mn]=25 [Fe]=26
        [Co]=27 [Ni]=28 [Cu]=29 [Zn]=30 [Ga]=31 [Ge]=32 [As]=33 [Se]=34
        [Br]=35 [Kr]=36 [Rb]=37 [Sr]=38 [Y]=39 [Zr]=40 [Nb]=41 [Mo]=42
        [Tc]=43 [Ru]=44 [Rh]=45 [Pd]=46 [Ag]=47 [Cd]=48 [In]=49 [Sn]=50
        [Sb]=51 [Te]=52 [I]=53 [Xe]=54 [Cs]=55 [Ba]=56 [La]=57 [Ce]=58
        [Pr]=59 [Nd]=60 [Pm]=61 [Sm]=62 [Eu]=63 [Gd]=64 [Tb]=65 [Dy]=66
        [Ho]=67 [Er]=68 [Tm]=69 [Yb]=70 [Lu]=71 [Hf]=72 [Ta]=73 [W]=74
        [Re]=75 [Os]=76 [Ir]=77 [Pt]=78 [Au]=79 [Hg]=80 [Tl]=81 [Pb]=82
        [Bi]=83 [Po]=84 [At]=85 [Rn]=86 [Fr]=87 [Ra]=88 [Ac]=89 [Th]=90
        [Pa]=91 [U]=92 [Np]=93 [Pu]=94
    )

    local total=0
    local line_num=0

    while IFS= read -r line; do
        ((line_num++))
        [ $line_num -le 2 ] && continue
        local element=$(echo "$line" | awk '{print $1}')
        [ -z "$element" ] && continue
        element=$(echo "$element" | sed 's/\(.\)\(.*\)/\U\1\L\2/')
        [ -n "${ATOMIC_NUMBERS[$element]}" ] && total=$((total + ${ATOMIC_NUMBERS[$element]}))
    done < "$xyz_file"

    echo $total
}

# =============================================================================
# Main script
# =============================================================================
echo -e "${BG_GREEN}${WHITE}====================================${NC}"
echo -e "${BG_GREEN}${WHITE}     ORCA 6.1 Workflow Setup        ${NC}"
echo -e "${BG_GREEN}${WHITE}====================================${NC}"
echo -e "ORCA Path: ${CYAN}${ORCA_PATH}${NC}"
echo ""

# Check SLURM cluster resources and select number of cores
check_slurm_resources
echo ""

# Check coord.xyz
if [ ! -f "coord.xyz" ]; then
    echo -e "${BG_RED}${WHITE} ERROR: coord.xyz not found in current directory. ${NC}"
    exit 1
fi
echo -e "${GREEN}[OK]${NC} coord.xyz found"

# =============================================================================
# Charge input
# =============================================================================
echo -e "${YELLOW}====================================${NC}"
read -p "Charge [default: 0]: " CHARGE_INPUT
CHARGE=${CHARGE_INPUT:-0}
echo -e "Charge: ${CYAN}${CHARGE}${NC}"

# =============================================================================
# Multiplicity input (auto estimation)
# =============================================================================
read -p "Multiplicity [Enter=auto]: " SPIN_INPUT

if [ -z "$SPIN_INPUT" ]; then
    echo "Estimating multiplicity..."
    TOTAL_ELECTRONS=$(get_total_electrons "coord.xyz")
    ELECTRONS_WITH_CHARGE=$((TOTAL_ELECTRONS - CHARGE))

    if [ $((ELECTRONS_WITH_CHARGE % 2)) -eq 0 ]; then
        SPIN=1
    else
        SPIN=2
    fi
    echo -e "Total electrons: ${TOTAL_ELECTRONS}, after charge: ${ELECTRONS_WITH_CHARGE}"
    echo -e "Estimated multiplicity: ${CYAN}${SPIN}${NC}"
else
    SPIN=$SPIN_INPUT
    echo -e "Multiplicity: ${CYAN}${SPIN}${NC}"
fi



# =============================================================================
# Calculation type
# =============================================================================
echo -e "${YELLOW}====================================${NC}"
echo "Calculation type:"
echo "  1) Opt > Freq"
echo "  2) Opt > Freq > SP (Compound)"
echo "  3) TS Opt > Freq"
echo "  4) TS Opt > Freq > SP (Compound)"
echo "  5) IRC"
echo "  6) Single Point"
echo "  7) Freq only"
echo "  8) SP Correction only (def2-TZVPPD)"

USE_COMPOUND=false
while true; do
    read -p "Select [1-8]: " CALC_CHOICE
    case $CALC_CHOICE in
        1) CALC_TYPE="OPT Freq"; CALC_NAME="opt_freq"; break ;;
        2) CALC_TYPE="OPT Freq"; CALC_NAME="opt_freq_sp"; USE_COMPOUND=true; break ;;
        3) CALC_TYPE="OptTS Freq"; CALC_NAME="tsopt_freq"; break ;;
        4) CALC_TYPE="OptTS Freq"; CALC_NAME="tsopt_freq_sp"; USE_COMPOUND=true; break ;;
        5) CALC_TYPE="IRC"; CALC_NAME="irc"; break ;;
        6) CALC_TYPE=""; CALC_NAME="sp"; break ;;
        7) CALC_TYPE="Freq"; CALC_NAME="freq"; break ;;
        8) CALC_TYPE=""; CALC_NAME="sp_correction"; SP_CORRECTION=true; break ;;
        *) echo -e "${BG_RED}${WHITE} Invalid choice. ${NC}" ;;
    esac
done
echo -e "Selected: ${CYAN}${CALC_NAME}${NC}"
[ "$USE_COMPOUND" = true ] && echo -e "${BG_MAGENTA}${WHITE} Compound Script: Freq(def2-SVP) > SP(def2-TZVPPD) ${NC}"

# =============================================================================
# Functional selection
# =============================================================================
echo -e "${YELLOW}====================================${NC}"
echo "Functional:"
echo "  1) r2SCAN-3c (fast & accurate)"
echo "  2) B3LYP-D3BJ [default]"
echo "  3) r2SCAN + D4"

while true; do
    read -p "Select [1-3, default=2]: " FUNC_CHOICE
    FUNC_CHOICE=${FUNC_CHOICE:-2}
    case $FUNC_CHOICE in
        1)
            FUNCTIONAL="r2SCAN-3c"
            DISPERSION=""
            USE_RI=""
            BASIS=""
            break
            ;;
        2)
            FUNCTIONAL="B3LYP"
            DISPERSION="D3BJ"
            USE_RI="def2/J RIJCOSX"
            BASIS="def2-SVP"
            break
            ;;
        3)
            FUNCTIONAL="r2SCAN"
            DISPERSION="D4"
            USE_RI="def2/J RIJCOSX"
            BASIS="def2-SVP"
            break
            ;;
        *) echo -e "${BG_RED}${WHITE} Invalid choice. ${NC}" ;;
    esac
done

# SP Correction: change basis to def2-TZVPPD
if [ "$SP_CORRECTION" = true ]; then
    if [ -n "$BASIS" ]; then
        BASIS="def2-TZVPPD"
    fi
    echo -e "${BG_MAGENTA}${WHITE} SP Correction mode: basis set = def2-TZVPPD ${NC}"
fi

echo -e "Functional: ${CYAN}${FUNCTIONAL}${NC}"
[ -n "$DISPERSION" ] && echo -e "Dispersion: ${CYAN}${DISPERSION}${NC}"
[ -n "$BASIS" ] && echo -e "Basis: ${CYAN}${BASIS}${NC}"

# =============================================================================
# Solvent selection
# =============================================================================
echo -e "${YELLOW}====================================${NC}"
echo "Solvent:"
echo "  1) DMF"
echo "  2) THF"
echo "  3) Water"
echo "  4) Acetonitrile"
echo "  5) None (gas phase)"

DEFAULT_SOL_NUM=1
[ "$DEFAULT_SOLVENT" = "THF" ] && DEFAULT_SOL_NUM=2
[ "$DEFAULT_SOLVENT" = "Water" ] && DEFAULT_SOL_NUM=3
[ "$DEFAULT_SOLVENT" = "Acetonitrile" ] && DEFAULT_SOL_NUM=4
[ -z "$DEFAULT_SOLVENT" ] && DEFAULT_SOL_NUM=5

while true; do
    read -p "Select [1-5, default=${DEFAULT_SOL_NUM}]: " SOLVENT_CHOICE
    SOLVENT_CHOICE=${SOLVENT_CHOICE:-$DEFAULT_SOL_NUM}
    case $SOLVENT_CHOICE in
        1) SOLVENT="CPCM(DMF)"; break ;;
        2) SOLVENT="CPCM(THF)"; break ;;
        3) SOLVENT="CPCM(Water)"; break ;;
        4) SOLVENT="CPCM(Acetonitrile)"; break ;;
        5) SOLVENT=""; break ;;
        *) echo -e "${BG_RED}${WHITE} Invalid choice. ${NC}" ;;
    esac
done
[ -n "$SOLVENT" ] && echo -e "Solvent: ${CYAN}${SOLVENT}${NC}"

# =============================================================================
# MO saving (ORCA 6 syntax)
# =============================================================================
echo -e "${YELLOW}====================================${NC}"
read -p "Save MO? [y/n, default=n]: " SAVE_MO

OUTPUT_BLOCK=""
if [[ "$SAVE_MO" =~ ^[Yy]$ ]]; then
    read -p "HOMO-N to LUMO+N, N= [default: ${DEFAULT_MO_RANGE}]: " N_ORBITALS
    N_ORBITALS=${N_ORBITALS:-$DEFAULT_MO_RANGE}

    OUTPUT_BLOCK="%output
  Print[P_MOs] 1
  PlotOrbs[HOMO-${N_ORBITALS}:LUMO+${N_ORBITALS}]
  PlotFormat Cube
  PlotNGridX ${CUBE_GRID_X}
  PlotNGridY ${CUBE_GRID_Y}
  PlotNGridZ ${CUBE_GRID_Z}
end
"
    echo -e "MO: ${CYAN}HOMO-${N_ORBITALS} ~ LUMO+${N_ORBITALS}${NC}, Grid: ${CUBE_GRID_X}x${CUBE_GRID_Y}x${CUBE_GRID_Z}"
fi

# =============================================================================
# Build keywords
# =============================================================================
KEYWORDS="! UKS ${FUNCTIONAL}"
[ -n "$BASIS" ] && KEYWORDS="${KEYWORDS} ${BASIS}"
[ -n "$USE_RI" ] && KEYWORDS="${KEYWORDS} ${USE_RI}"
[ -n "$DISPERSION" ] && KEYWORDS="${KEYWORDS} ${DISPERSION}"
[ -n "$CALC_TYPE" ] && KEYWORDS="${KEYWORDS} ${CALC_TYPE}"
[ -n "$SOLVENT" ] && KEYWORDS="${KEYWORDS} ${SOLVENT}"
KEYWORDS="${KEYWORDS} TightOpt TightSCF"

if [ "$FUNCTIONAL" != "r2SCAN-3c" ]; then
    # Use fine grid for SP correction, default grid for others
    if [ "$SP_CORRECTION" = true ]; then
        KEYWORDS="${KEYWORDS} ${DFT_GRID_SP}"
    else
        KEYWORDS="${KEYWORDS} ${DFT_GRID}"
    fi
fi

# =============================================================================
# %geom block (for TS optimization)
# =============================================================================
GEOM_BLOCK=""
if [[ "$CALC_TYPE" == *"OptTS"* ]]; then
    echo -e "${YELLOW}====================================${NC}"
    read -p "TS Mode index [default: 0]: " TS_MODE
    TS_MODE=${TS_MODE:-0}

    GEOM_BLOCK="%geom
  Calc_Hess true
  Recalc_Hess 3
  Trust 0.2
  MaxIter 300
  TS_Mode {M ${TS_MODE}}
end
"
fi

# =============================================================================
# IRC block
# =============================================================================
IRC_BLOCK=""
if [[ "$CALC_TYPE" == "IRC" ]]; then
    IRC_BLOCK="%irc
  MaxIter 100
  Direction Both
end
"
    echo -e "${BG_MAGENTA}${WHITE} WARNING: IRC requires .hess file from previous TS calculation. ${NC}"
fi

# =============================================================================
# Generate filename (includes current folder name)
# =============================================================================
FOLDER_NAME="${PWD##*/}"
if [ "$SP_CORRECTION" = true ]; then
    OUTPUT_NAME="${FOLDER_NAME}_${CALC_NAME}_${FUNCTIONAL}_def2-TZVPPD"
else
    OUTPUT_NAME="${FOLDER_NAME}_${CALC_NAME}_${FUNCTIONAL}"
    [ -n "$BASIS" ] && OUTPUT_NAME="${OUTPUT_NAME}_${BASIS}"
fi
OUTPUT_NAME=$(echo "$OUTPUT_NAME" | tr ' ' '_' | tr -cd '[:alnum:]_-')

# =============================================================================
# Generate input.inp
# =============================================================================
echo -e "${YELLOW}====================================${NC}"
echo -e "Creating: ${GREEN}${OUTPUT_NAME}.inp${NC}"

if [ "$USE_COMPOUND" = true ]; then
    # Compound Script: Step 1 (Opt+Freq with small basis) > Step 2 (SP with large basis)

    # Build Step 1 keywords (Opt + Freq with small basis)
    STEP1_KEYWORDS="! UKS ${FUNCTIONAL}"
    if [ "$FUNCTIONAL" = "r2SCAN-3c" ]; then
        # r2SCAN-3c has built-in basis
        STEP1_KEYWORDS="${STEP1_KEYWORDS}"
    else
        STEP1_KEYWORDS="${STEP1_KEYWORDS} def2-SVP"
        [ -n "$USE_RI" ] && STEP1_KEYWORDS="${STEP1_KEYWORDS} ${USE_RI}"
    fi
    [ -n "$DISPERSION" ] && STEP1_KEYWORDS="${STEP1_KEYWORDS} ${DISPERSION}"
    STEP1_KEYWORDS="${STEP1_KEYWORDS} ${CALC_TYPE}"
    [ -n "$SOLVENT" ] && STEP1_KEYWORDS="${STEP1_KEYWORDS} ${SOLVENT}"
    STEP1_KEYWORDS="${STEP1_KEYWORDS} TightOpt TightSCF"
    [ "$FUNCTIONAL" != "r2SCAN-3c" ] && STEP1_KEYWORDS="${STEP1_KEYWORDS} ${DFT_GRID}"

    # Build Step 2 keywords (SP with large basis and fine grid)
    STEP2_KEYWORDS="! UKS ${FUNCTIONAL}"
    if [ "$FUNCTIONAL" = "r2SCAN-3c" ]; then
        # For r2SCAN-3c, switch to r2SCAN with explicit large basis for SP
        STEP2_KEYWORDS="! UKS r2SCAN def2-TZVPPD def2/J RIJCOSX D4"
    else
        STEP2_KEYWORDS="${STEP2_KEYWORDS} def2-TZVPPD"
        [ -n "$USE_RI" ] && STEP2_KEYWORDS="${STEP2_KEYWORDS} ${USE_RI}"
        [ -n "$DISPERSION" ] && STEP2_KEYWORDS="${STEP2_KEYWORDS} ${DISPERSION}"
    fi
    [ -n "$SOLVENT" ] && STEP2_KEYWORDS="${STEP2_KEYWORDS} ${SOLVENT}"
    STEP2_KEYWORDS="${STEP2_KEYWORDS} TightSCF ${DFT_GRID_SP}"

    cat > "${OUTPUT_NAME}.inp" << EOF
# ORCA 6.1 Compound Script
# Generated: $(date)
# Functional: ${FUNCTIONAL}
# Step 1: ${CALC_TYPE} with def2-SVP (or r2SCAN-3c basis)
# Step 2: SP with def2-TZVPPD

# =============================================================================
# Step 1: Optimization + Frequency (small basis)
# =============================================================================
${STEP1_KEYWORDS}

${GEOM_BLOCK}${OUTPUT_BLOCK}%pal
  nprocs ${NPROCS}
end

%maxcore ${MAXCORE}

*xyzfile ${CHARGE} ${SPIN} coord.xyz

# =============================================================================
# Step 2: Single Point Correction (large basis - def2-TZVPPD)
# =============================================================================
\$new_job
${STEP2_KEYWORDS}

%pal
  nprocs ${NPROCS}
end

%maxcore ${MAXCORE}

*xyzfile ${CHARGE} ${SPIN} ${OUTPUT_NAME}.xyz
EOF

    echo -e "${GREEN}[OK]${NC} ${OUTPUT_NAME}.inp created (Compound Script)"
    echo -e "${CYAN}Step 1:${NC} ${CALC_TYPE} with def2-SVP"
    echo -e "${CYAN}Step 2:${NC} SP with def2-TZVPPD"

else
    # Standard single-step calculation
    cat > "${OUTPUT_NAME}.inp" << EOF
# ORCA 6.1 Input File
# Generated: $(date)
# Functional: ${FUNCTIONAL}
# Basis: ${BASIS:-"included in method"}
# Calculation: ${CALC_NAME}

${KEYWORDS}

${GEOM_BLOCK}${IRC_BLOCK}${OUTPUT_BLOCK}%pal
  nprocs ${NPROCS}
end

%maxcore ${MAXCORE}

*xyzfile ${CHARGE} ${SPIN} coord.xyz
EOF

    echo -e "${GREEN}[OK]${NC} ${OUTPUT_NAME}.inp created"
fi

# =============================================================================
# Generate SLURM submission script
# =============================================================================
PARTITION_LINE=""
[ -n "$SLURM_PARTITION" ] && PARTITION_LINE="#SBATCH --partition=${SLURM_PARTITION}"

# Determine if energy extraction is needed
EXTRACT_ENERGY=false
if [ "$USE_COMPOUND" = true ] || [[ "$CALC_TYPE" == *"Freq"* ]]; then
    EXTRACT_ENERGY=true
fi

cat > "submit_${OUTPUT_NAME}.sh" << EOF
#!/bin/bash
#SBATCH --job-name=${OUTPUT_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks=${NPROCS}
#SBATCH --cpus-per-task=1
#SBATCH --mem=${SLURM_MEM}
#SBATCH --time=${SLURM_TIME}
${PARTITION_LINE}
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

export PATH="${ORCA_PATH}:\${PATH}"
export LD_LIBRARY_PATH="${ORCA_PATH}:\${LD_LIBRARY_PATH}"

cd \${SLURM_SUBMIT_DIR}

echo "========================================"
echo "ORCA 6.1 Calculation"
echo "========================================"
echo "Job:     \${SLURM_JOB_NAME}"
echo "ID:      \${SLURM_JOB_ID}"
echo "Node:    \$(hostname)"
echo "Start:   \$(date)"
echo "========================================"

${ORCA_PATH}/orca ${OUTPUT_NAME}.inp > ${OUTPUT_NAME}.out 2>&1

echo "========================================"
echo "End:     \$(date)"
echo "========================================"
EOF

# Add energy extraction for Freq and Compound calculations
if [ "$EXTRACT_ENERGY" = true ]; then
    cat >> "submit_${OUTPUT_NAME}.sh" << 'ENERGY_SCRIPT'

# =============================================================================
# Energy Extraction and Summary
# =============================================================================
OUT_FILE="${OUTPUT_NAME}.out"
SUMMARY_FILE="${OUTPUT_NAME}_energy_summary.txt"

echo ""
echo "========================================"
echo "Extracting energies..."
echo "========================================"

# Extract energies - search from bottom of file and show full lines
E_ELECT_LINE=$(tac "$OUT_FILE" | grep -m1 "Electronic energy")
TOTAL_ENERGY_LINE=$(tac "$OUT_FILE" | grep -m1 "Total energy")
DISP_CORR_LINE=$(tac "$OUT_FILE" | grep -m1 "Dispersion correction")
ZPE_LINE=$(tac "$OUT_FILE" | grep -m1 "Non-thermal (ZPE) correction")
THERMAL_CORR_LINE=$(tac "$OUT_FILE" | grep -m1 "Total thermal correction")
ENTROPY_LINE=$(tac "$OUT_FILE" | grep -m1 "Final entropy term")

# Extract numeric values for calculations
E_ELECT=$(echo "$E_ELECT_LINE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}' | head -1)
TOTAL_ENERGY=$(echo "$TOTAL_ENERGY_LINE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}' | head -1)
DISP_CORR=$(echo "$DISP_CORR_LINE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}' | head -1)
ZPE=$(echo "$ZPE_LINE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}' | head -1)
G_CORR=$(echo "$THERMAL_CORR_LINE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}' | head -1)
ENTROPY=$(echo "$ENTROPY_LINE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}' | head -1)

# Also extract FINAL SINGLE POINT ENERGY for backward compatibility
FINAL_SP=$(grep "FINAL SINGLE POINT ENERGY" "$OUT_FILE" | tail -1 | awk '{print $5}')

# Extract Gibbs free energy from thermochemistry
G_THERMO=$(tac "$OUT_FILE" | grep -m1 "Final Gibbs free energy" | awk '{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}' | head -1)



ENERGY_SCRIPT

    # Add compound-specific extraction
    if [ "$USE_COMPOUND" = true ]; then
        cat >> "submit_${OUTPUT_NAME}.sh" << 'COMPOUND_SCRIPT'

# For Compound Script: Extract energies from both steps

# Find boundary between Step 1 and Step 2
STEP2_START=$(grep -n "\$new_job" "$OUT_FILE" | head -1 | cut -d: -f1)

if [ -z "$STEP2_START" ]; then
    # Fallback: look for "Step 2:" marker or second calculation start
    STEP2_START=$(grep -n "Step 2:" "$OUT_FILE" | head -1 | cut -d: -f1)
fi

# Extract Step 1 energies (from beginning to Step 2 start)
E_ELECT_STEP1_LINE=$(sed -n "1,${STEP2_START}p" "$OUT_FILE" | tac | grep -m1 "Electronic energy")
ZPE_STEP1_LINE=$(sed -n "1,${STEP2_START}p" "$OUT_FILE" | tac | grep -m1 "Non-thermal (ZPE) correction")
G_CORR_STEP1_LINE=$(sed -n "1,${STEP2_START}p" "$OUT_FILE" | tac | grep -m1 "Total thermal correction")
G_THERMO_STEP1_LINE=$(sed -n "1,${STEP2_START}p" "$OUT_FILE" | tac | grep -m1 "Final Gibbs free energy")

E_FREQ=$(echo "$E_ELECT_STEP1_LINE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}' | head -1)
ZPE_STEP1=$(echo "$ZPE_STEP1_LINE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}' | head -1)
G_CORR_STEP1=$(echo "$G_CORR_STEP1_LINE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}' | head -1)
G_THERMO_STEP1=$(echo "$G_THERMO_STEP1_LINE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}' | head -1)

# Extract Step 2 energies (from Step 2 start to end)
E_ELECT_STEP2_LINE=$(sed -n "${STEP2_START},\$p" "$OUT_FILE" | tac | grep -m1 "Electronic energy")
E_SP=$(echo "$E_ELECT_STEP2_LINE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}' | head -1)

# Calculate corrected Gibbs energy: G = E(SP) + G_corr
if [ -n "$E_SP" ] && [ -n "$G_CORR_STEP1" ]; then
    G_CORRECTED=$(echo "$E_SP + $G_CORR_STEP1" | bc -l)
fi


COMPOUND_SCRIPT
    fi

    # Add summary output
    cat >> "submit_${OUTPUT_NAME}.sh" << 'SUMMARY_SCRIPT'
# Write summary file
{
    echo "========================================"
    echo "     ORCA Energy Summary"
    echo "========================================"
    echo "Output file: $OUT_FILE"
    echo "Generated:   $(date)"
    echo "========================================"
    echo ""
SUMMARY_SCRIPT

    if [ "$USE_COMPOUND" = true ]; then
        cat >> "submit_${OUTPUT_NAME}.sh" << 'COMPOUND_SUMMARY'
    echo "[Compound Script Results]"
    echo ""
    echo "Step 1 (Opt+Freq / def2-SVP):"
    if [ -n "$E_ELECT_STEP1_LINE" ]; then
        echo "  $E_ELECT_STEP1_LINE"
    fi
    if [ -n "$ZPE_STEP1_LINE" ]; then
        echo "  $ZPE_STEP1_LINE"
    fi
    if [ -n "$G_CORR_STEP1_LINE" ]; then
        echo "  $G_CORR_STEP1_LINE"
    fi
    if [ -n "$G_THERMO_STEP1_LINE" ]; then
        echo "  $G_THERMO_STEP1_LINE"
    fi
    echo ""
    echo "Step 2 (SP / def2-TZVPPD):"
    [ -n "$E_ELECT_STEP2_LINE" ] && echo "  $E_ELECT_STEP2_LINE"
    echo ""
    echo "----------------------------------------"
    echo "[Final Corrected Gibbs Energy]"
    echo "  G = E(SP/TZVPPD) + G_corr(Freq/SVP)"
    [ -n "$G_CORRECTED" ] && echo "  G = $G_CORRECTED Eh"
    echo "----------------------------------------"
COMPOUND_SUMMARY
    else
        cat >> "submit_${OUTPUT_NAME}.sh" << 'SINGLE_SUMMARY'

SINGLE_SUMMARY
    fi
    echo "[Single Calculation Results]"
    echo ""
    if [ -n "$E_ELECT_LINE" ]; then
        echo "  $E_ELECT_LINE"
    fi
    if [ -n "$TOTAL_ENERGY_LINE" ]; then
        echo "  $TOTAL_ENERGY_LINE"
    fi
    if [ -n "$DISP_CORR_LINE" ]; then
        echo "  $DISP_CORR_LINE"
    fi
    if [ -n "$ZPE_LINE" ]; then
        echo "  $ZPE_LINE"
    fi
    if [ -n "$THERMAL_CORR_LINE" ]; then
        echo "  $THERMAL_CORR_LINE"
    fi
    if [ -n "$ENTROPY_LINE" ]; then
        echo "  $ENTROPY_LINE"
    fi
    if [ -n "$G_THERMO" ]; then
        echo ""
        echo "  Final Gibbs free energy:              $G_THERMO Eh"
    fi
    cat >> "submit_${OUTPUT_NAME}.sh" << 'END_SUMMARY'
    echo ""
    echo "========================================"
} > "$SUMMARY_FILE"

# Display summary
cat "$SUMMARY_FILE"
echo ""
echo "Summary saved to: $SUMMARY_FILE"
END_SUMMARY

    # Replace placeholder with actual output name
    sed -i "s/\${OUTPUT_NAME}/${OUTPUT_NAME}/g" "submit_${OUTPUT_NAME}.sh"
fi

chmod +x "submit_${OUTPUT_NAME}.sh"
echo -e "${GREEN}[OK]${NC} submit_${OUTPUT_NAME}.sh created"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo -e "${BG_GREEN}${WHITE}====================================${NC}"
echo -e "${BG_GREEN}${WHITE}            Summary                 ${NC}"
echo -e "${BG_GREEN}${WHITE}====================================${NC}"
echo -e "Charge:     ${CYAN}${CHARGE}${NC}"
echo -e "Mult:       ${CYAN}${SPIN}${NC}"
echo -e "Calc:       ${CYAN}${CALC_NAME}${NC}"
echo -e "Functional: ${CYAN}${FUNCTIONAL}${NC}"
[ -n "$DISPERSION" ] && echo -e "Dispersion: ${CYAN}${DISPERSION}${NC}"
if [ "$USE_COMPOUND" = true ]; then
    echo -e "Basis:      ${CYAN}def2-SVP (Opt/Freq) > def2-TZVPPD (SP)${NC}"
else
    [ -n "$BASIS" ] && echo -e "Basis:      ${CYAN}${BASIS}${NC}"
fi
[ -n "$SOLVENT" ] && echo -e "Solvent:    ${CYAN}${SOLVENT}${NC}"
echo -e "Cores:      ${CYAN}${NPROCS}${NC}"
echo -e "Memory:     ${CYAN}${MAXCORE} MB/core${NC}"
if [ "$USE_COMPOUND" = true ]; then
    echo ""
    echo -e "${BG_MAGENTA}${WHITE} Compound Script Workflow: ${NC}"
    echo -e "  Step 1: ${CALC_TYPE} (def2-SVP)"
    echo -e "  Step 2: SP Correction (def2-TZVPPD)"
    echo -e "  G = E(SP/TZVPPD) + G_corr(Freq/SVP)"
fi
echo ""
echo -e "${BG_GREEN}${WHITE}====================================${NC}"
echo -e "Run: ${CYAN}sbatch submit_${OUTPUT_NAME}.sh${NC}"
echo -e "${BG_GREEN}${WHITE}====================================${NC}"

# =============================================================================
# Auto submit option
# =============================================================================
read -p "Submit now? [y/n, default=y]: " SUBMIT_NOW
SUBMIT_NOW=${SUBMIT_NOW:-y}
if [[ "$SUBMIT_NOW" =~ ^[Yy]$ ]]; then
    sbatch "submit_${OUTPUT_NAME}.sh"

fi