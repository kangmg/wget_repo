#!/bin/bash
# =============================================================================
#                    ORCA 6.1 Workflow Script
# =============================================================================
# SLURM-based / 6-core environment / Auto charge/multiplicity estimation
# =============================================================================

# =============================================================================
#                         USER CONFIGURATION
# =============================================================================
# Modify these variables according to your environment

# ORCA installation path (REQUIRED)
ORCA_PATH="/path/to/orca_6.1"

# Number of CPU cores per job
NPROCS=6

# Memory per core (MB) - Total memory = NPROCS * MAXCORE
MAXCORE=5000

# Default solvent for CPCM (DMF, THF, Water, Acetonitrile, or empty for gas phase)
DEFAULT_SOLVENT="DMF"

# MO Cube file settings
CUBE_GRID_X=80
CUBE_GRID_Y=80
CUBE_GRID_Z=80
DEFAULT_MO_RANGE=3    # HOMO-N to LUMO+N

# SLURM settings
SLURM_PARTITION=""           # Leave empty for default partition
SLURM_TIME="72:00:00"        # Max walltime
SLURM_MEM="32G"              # Total memory request

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
echo -e "Cores: ${CYAN}${NPROCS}${NC}, Memory: ${CYAN}$((NPROCS * MAXCORE / 1000))GB${NC}"
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
echo "  1) Opt + Freq"
echo "  2) TS Opt + Freq"
echo "  3) IRC"
echo "  4) Single Point"
echo "  5) Freq only"
echo "  6) SP Correction (def2-TZVPPD)"

while true; do
    read -p "Select [1-6]: " CALC_CHOICE
    case $CALC_CHOICE in
        1) CALC_TYPE="OPT Freq"; CALC_NAME="opt_freq"; break ;;
        2) CALC_TYPE="OptTS Freq"; CALC_NAME="tsopt_freq"; break ;;
        3) CALC_TYPE="IRC"; CALC_NAME="irc"; break ;;
        4) CALC_TYPE=""; CALC_NAME="sp"; break ;;
        5) CALC_TYPE="Freq"; CALC_NAME="freq"; break ;;
        6) CALC_TYPE=""; CALC_NAME="sp_correction"; SP_CORRECTION=true; break ;;
        *) echo -e "${BG_RED}${WHITE} Invalid choice. ${NC}" ;;
    esac
done
echo -e "Selected: ${CYAN}${CALC_NAME}${NC}"

# =============================================================================
# Functional selection
# =============================================================================
echo -e "${YELLOW}====================================${NC}"
echo "Functional:"
echo "  1) r2SCAN-3c (default - fast & accurate)"
echo "  2) B3LYP-D3BJ"
echo "  3) B3LYP (no dispersion)"
echo "  4) r2SCAN + D4"

while true; do
    read -p "Select [1-4, default=1]: " FUNC_CHOICE
    FUNC_CHOICE=${FUNC_CHOICE:-1}
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
            FUNCTIONAL="B3LYP"
            DISPERSION=""
            USE_RI="def2/J RIJCOSX"
            BASIS="def2-SVP"
            break
            ;;
        4)
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
read -p "Save MO? [y/N]: " SAVE_MO

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
    KEYWORDS="${KEYWORDS} Grid5 FinalGrid6"
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
# Generate filename
# =============================================================================
if [ "$SP_CORRECTION" = true ]; then
    OUTPUT_NAME="${CALC_NAME}_${FUNCTIONAL}_def2-TZVPPD"
else
    OUTPUT_NAME="${CALC_NAME}_${FUNCTIONAL}"
    [ -n "$BASIS" ] && OUTPUT_NAME="${OUTPUT_NAME}_${BASIS}"
fi
OUTPUT_NAME=$(echo "$OUTPUT_NAME" | tr ' ' '_' | tr -cd '[:alnum:]_-')

# =============================================================================
# Generate input.inp
# =============================================================================
echo -e "${YELLOW}====================================${NC}"
echo -e "Creating: ${GREEN}${OUTPUT_NAME}.inp${NC}"

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

# =============================================================================
# Generate SLURM submission script
# =============================================================================
PARTITION_LINE=""
[ -n "$SLURM_PARTITION" ] && PARTITION_LINE="#SBATCH --partition=${SLURM_PARTITION}"

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
[ -n "$BASIS" ] && echo -e "Basis:      ${CYAN}${BASIS}${NC}"
[ -n "$SOLVENT" ] && echo -e "Solvent:    ${CYAN}${SOLVENT}${NC}"
echo -e "Cores:      ${CYAN}${NPROCS}${NC}"
echo -e "Memory:     ${CYAN}${MAXCORE} MB/core${NC}"
echo ""
echo -e "${BG_GREEN}${WHITE}====================================${NC}"
echo -e "Run: ${CYAN}sbatch submit_${OUTPUT_NAME}.sh${NC}"
echo -e "${BG_GREEN}${WHITE}====================================${NC}"

# =============================================================================
# Auto submit option
# =============================================================================
read -p "Submit now? [y/N]: " SUBMIT_NOW
if [[ "$SUBMIT_NOW" =~ ^[Yy]$ ]]; then
    sbatch "submit_${OUTPUT_NAME}.sh"
fi