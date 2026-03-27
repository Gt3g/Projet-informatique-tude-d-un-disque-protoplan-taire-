# ================================================================
#  Makefile — Simulation de formation d'un disque protoplanétaire
#
#  Cibles disponibles :
#    make          → compile l'exécutable (mode release)
#    make run      → compile puis lance la simulation
#    make debug    → compile avec symboles de debug + sanitizers
#    make valgrind → compile en debug puis lance sous Valgrind
#    make plot     → lance le script de visualisation Python
#    make clean    → supprime les .o, l'exécutable et les fichiers de sortie
#    make help     → affiche ce résumé
#
#  Prérequis :
#    g++ ≥ 5 (C++14), make, python3 + numpy + matplotlib (pour make plot)
#
#  Tous les paramètres de simulation sont dans params.h.
#  Recompiler avec 'make' après toute modification de params.h.
# ================================================================

# ── Compilateur et options ────────────────────────────────────────
CXX      := g++
CXXFLAGS := -std=c++14 -Wall -Wextra -O2
LDFLAGS  :=

# ── Sources, objets, exécutable ───────────────────────────────────
SRCS := main_nuage_etoile.cpp \
        particule.cpp         \
        nuage.cpp             \
        nuage_solver.cpp      \
        forces.cpp

OBJS   := $(SRCS:.cpp=.o)
TARGET := simulation_disque

# Script Python de visualisation (optionnel)
PLOT_SCRIPT := visualisation.py

# ── Cible par défaut : compilation release ────────────────────────
.PHONY: all
all: $(TARGET)
	@echo ""
	@echo "  ✓  Compilation réussie → ./$(TARGET)"
	@echo "  →  Lancer la simulation  : make run"
	@echo "  →  Visualiser les résultats : make plot"
	@echo ""

# ── Édition de liens ──────────────────────────────────────────────
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# ── Règle générique .cpp → .o ─────────────────────────────────────
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# ── Dépendances explicites des en-têtes ──────────────────────────
#    (garantit la recompilation si un .h est modifié)
main_nuage_etoile.o : main_nuage_etoile.cpp params.h nuage.h nuage_solver.h
particule.o         : particule.cpp particule.h
nuage.o             : nuage.cpp nuage.h particule.h
nuage_solver.o      : nuage_solver.cpp nuage_solver.h nuage.h forces.h params.h
forces.o            : forces.cpp forces.h nuage.h nuage_solver.h params.h

# ── Lancer la simulation ──────────────────────────────────────────
.PHONY: run
run: $(TARGET)
	@echo "  →  Démarrage de la simulation..."
	./$(TARGET)

# ── Mode debug (avec AddressSanitizer + UBSan) ────────────────────
#    Détecte : débordements mémoire, comportements indéfinis, etc.
#    Usage : make debug && ./simulation_disque
.PHONY: debug
debug: CXXFLAGS := -std=c++14 -Wall -Wextra -g -O0 \
                   -fsanitize=address,undefined \
                   -fno-omit-frame-pointer
debug: LDFLAGS  := -fsanitize=address,undefined
debug: $(TARGET)
	@echo ""
	@echo "  ✓  Build debug → ./$(TARGET)  (sanitizers activés)"
	@echo ""

# ── Valgrind (mémoire) ────────────────────────────────────────────
#    Nécessite : apt install valgrind
#    Plus lent qu'ASan mais compatible sans recompilation spéciale.
.PHONY: valgrind
valgrind: CXXFLAGS := -std=c++14 -Wall -Wextra -g -O0
valgrind: LDFLAGS  :=
valgrind: $(TARGET)
	valgrind --leak-check=full --track-origins=yes \
	         --error-exitcode=1 ./$(TARGET)

# ── Visualisation Python ──────────────────────────────────────────
#    Lance visualisation.py si présent ; affiche un message sinon.
.PHONY: plot
plot:
	@if [ -f "$(PLOT_SCRIPT)" ]; then \
	    echo "  →  Lancement de $(PLOT_SCRIPT)..."; \
	    python3 $(PLOT_SCRIPT); \
	else \
	    echo "  ✗  $(PLOT_SCRIPT) introuvable."; \
	    echo "     Créez ce script ou adaptez le chemin dans le Makefile (variable PLOT_SCRIPT)."; \
	fi

# ── Nettoyage ─────────────────────────────────────────────────────
.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET) evolution_nuage.fich diagnostics.fich
	@echo "  ✓  Nettoyage terminé."

# ── Aide ──────────────────────────────────────────────────────────
.PHONY: help
help:
	@echo ""
	@echo "  Simulation de formation d'un disque protoplanétaire"
	@echo "  ─────────────────────────────────────────────────────"
	@echo "  make           Compile en mode release (-O2)"
	@echo "  make run       Compile + lance la simulation"
	@echo "  make debug     Compile avec AddressSanitizer/UBSan"
	@echo "  make valgrind  Lance la simulation sous Valgrind"
	@echo "  make plot      Lance le script Python de visualisation"
	@echo "  make clean     Supprime objets, exécutable et fichiers de sortie"
	@echo "  make help      Affiche cette aide"
	@echo ""
	@echo "  Paramètres : éditer params.h puis 'make'"
	@echo ""
