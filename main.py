import sys
from pathlib import Path
from typing import List, NamedTuple, Tuple

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QAction, QColor, QFont, QFontMetrics, QMouseEvent, QPainter
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMenuBar,
    QPushButton,
    QSpinBox,
    QSplitter,
    QStatusBar,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)


class SequenceBrowser(QWidget):
    """Widget central affichant la séquence nucléotidique avec zoom, pan et traçage."""

    position_changed = pyqtSignal(int, str)  # (position_1-based, base)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.sequence = ""

        # Configuration typographique
        font = QFont("Inconsolata", 11)
        font.setFixedPitch(True)
        self.setFont(font)
        self.setStyleSheet("background-color: #FAFAFA;")

        # Gestion du focus
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setFocus()

        # État de la vue
        self.scroll_x = 0
        self.zoom_factor = 1.0

    def load_fasta(self, filepath: str):
        """Lit un fichier FASTA et extrait la séquence brute."""
        if not Path(filepath).exists():
            print(f"Fichier introuvable : {filepath}")
            return

        with open(filepath, "r") as f:
            lines = [L.strip() for L in f.readlines()]

        # Extraction propre (ignore les headers et retours ligne)
        seq_lines = [line for line in lines if not line.startswith(">")]
        self.sequence = "".join(seq_lines).upper().replace(" ", "")

        # Ajouter focus
        self.scroll_x = 0
        self.setFocus()
        self.update()

    def _get_metrics(self):
        """Retourne la largeur d'une base et la hauteur de ligne en fonction du zoom."""
        fm = QFontMetrics(self.font())
        bw = max(fm.horizontalAdvance("A"), 18) * self.zoom_factor
        return bw, fm.height()

    def paintEvent(self, event):  # type: ignore[override]
        super().paintEvent(event)
        painter = QPainter(self)
        painter.setFont(self.font())

        bw, lh = self._get_metrics()
        step = bw + 2  # espacement entre bases

        # Zones visibles
        vp_width = self.width()
        start_idx = max(0, int(self.scroll_x / step))
        num_visible = int(vp_width / step) + 1
        end_idx = min(start_idx + num_visible, len(self.sequence))

        # Centrage vertical du texte
        y_pos = (self.height() - lh) // 2

        for i in range(start_idx, end_idx):
            x = (i * step) - self.scroll_x
            base = self.sequence[i]

            # Coloration ACGT
            colors = {
                "A": QColor("#E63946"),
                "T": QColor("#2A9D8F"),
                "C": QColor("#1D3557"),
                "G": QColor("#E9C46A"),
            }
            painter.setPen(colors.get(base, QColor("black")))

            # Correction pour centrage vertical avec QFontMetrics
            y_offset = y_pos + (lh // 2)
            painter.drawText(int(x), int(y_offset), base)

        painter.end()

    def mouseMoveEvent(self, event: QMouseEvent):  # type: ignore[override]
        """Calcule la position du curseur sur la séquence et émet un signal."""
        if not self.sequence:
            return

        x = event.position().x() + self.scroll_x
        step = self._get_metrics()[0] + 2
        idx = int(x / step)
        idx = max(0, min(idx, len(self.sequence) - 1))

        self.position_changed.emit(idx + 1, self.sequence[idx])  # Converti en 1-based
        super().mouseMoveEvent(event)

    def mousePressEvent(self, event):  # type: ignore[override]
        """Focus sur le widget sur clic"""
        self.setFocus()
        super().mousePressEvent(event)

    def wheelEvent(self, event):  # type: ignore[override]
        """Gestion du zoom molette (0.5x à 3.0x)."""
        delta = event.angleDelta().y()
        if delta > 0:
            self.zoom_factor = min(self.zoom_factor * 1.15, 3.0)
        else:
            self.zoom_factor = max(self.zoom_factor / 1.15, 0.5)
        self.update()
        event.accept()

    def keyPressEvent(self, event):  # type: ignore[override]
        """Navigation clavier (flèches gauche/droite)."""
        step_size = int(20 * self._get_metrics()[0])
        if event.key() == Qt.Key.Key_Left:
            self.scroll_x = max(0, self.scroll_x - step_size)
        elif event.key() == Qt.Key.Key_Right:
            self.scroll_x += step_size
        self.update()
        super().keyPressEvent(event)


class Feature(NamedTuple):
    """Représente une caractéristique génomique (CDS, RBS, Promoteur)."""

    type: str  # 'CDS', 'RBS', 'Promoter'
    start: int  # Position 1-based
    end: int  # Position 1-based
    strand: str  # '+' ou '-'
    sequence: str  # La séquence nucléotidique correspondante
    length: int  # Longueur en nucléotides


class GenomeAnalyzer:
    """Classe statique pour l'analyse de la séquence ADN."""

    STOP_CODONS = {"TAA", "TAG", "TGA"}

    @staticmethod
    def reverse_complement(seq: str) -> str:
        """Génère le brin complémentaire inversé (3'->5')."""
        trans = str.maketrans("ACGT", "TGCA")
        return seq[::-1].translate(trans)

    @staticmethod
    def find_orfs(
        sequence: str, min_length: int, forward: bool, reverse: bool
    ) -> List[Feature]:
        """
        Détecte les ORFs.

        Args:
            sequence: La séquence ADN brute.
            min_length: Longueur minimale en nucléotides.
            forward: Scanner le brin +.
            reverse: Scanner le brin -.

        Retourne:
            Liste des objets Feature correspondant aux CDS candidats.
        """
        orfs = []

        def scan_strand(strand_name: str, seq: str):
            for frame in range(3):
                i = frame
                while i < len(seq) - 2:
                    codon = seq[i : i + 3]

                    if codon == "ATG":
                        start_idx = i
                        end_idx = None

                        # Chercher le premier stop codon en aval à partir du prochain triplet
                        j = i + 3
                        while j < len(seq) - 2:
                            next_codon = seq[j : j + 3]
                            if next_codon in GenomeAnalyzer.STOP_CODONS:
                                end_idx = j + 2
                                break
                            j += 3

                        if end_idx is not None:
                            cds_length = end_idx - start_idx + 1
                            if cds_length >= min_length:
                                raw_seq = seq[start_idx : end_idx + 1]

                                # Si c'est le brin reverse, on doit retourner la séquence complémentaire

                                actual_seq = raw_seq

                                orfs.append(
                                    Feature(
                                        type="CDS",
                                        start=start_idx + 1,  # 1-based
                                        end=end_idx + 1,  # 1-based
                                        strand=strand_name,
                                        sequence=actual_seq,
                                        length=cds_length,
                                    )
                                )

                        i += 3
                    else:
                        i += 3

        # Brins à scanner
        if forward:
            scan_strand("+", sequence)

        if reverse:
            rev_seq = GenomeAnalyzer.reverse_complement(sequence)
            scan_strand("-", rev_seq)

        return orfs

    @staticmethod
    def find_motifs(
        sequence: str, motif: str, window_start: int, window_end: int
    ) -> List[Feature]:
        """
        Recherche un motif dans une fenêtre de position relative.

        Args:
            sequence: La séquence ADN (brin +).
            motif: Le motif à chercher (ex: TTGACA).
            window_start: Position absolue du début de la zone de recherche (1-based).
            window_end: Position absolue de la fin de la zone de recherche (1-based).

        Retourne:
            Liste des Features trouvés.
        """
        found = []
        search_seq = sequence[window_start - 1 : window_end]

        for i in range(len(search_seq) - len(motif) + 1):
            if search_seq[i : i + len(motif)] == motif:
                abs_pos = window_start + i  # Position absolue dans la séquence globale
                found.append(
                    Feature(
                        type="Promoter"
                        if "TTGACA" in motif or "TATAAT" in motif
                        else "RBS",
                        start=abs_pos,
                        end=abs_pos + len(motif) - 1,
                        strand="+",
                        sequence=motif,
                        length=len(motif),
                    )
                )
        return found

    @staticmethod
    def find_rbs_upstream(
        sequence: str,
        cds_feature: Feature,
        motif_pattern: str = "AGGAGGUAA",
        upstream_dist_min: int = 6,
        upstream_dist_max: int = 12,
    ) -> List[Feature]:
        """
        Cherche un site RBS (Shine-Dalgarno) en amont d'un CDS.

        Args:
            sequence: La séquence génomique complète (+).
            cds_feature: L'objet Feature du CDS trouvé.
            motif_pattern: Motif de recherche pour le RBS.

        Retourne:
            Liste des Features RBS trouvés.
        """
        rbs_list = []

        # Si le CDS est sur le brin -, le RBS est en amont du Start (donc vers les indices plus grands)
        # Si le CDS est sur le brin +, le RBS est en amont du Start (indices plus petits)

        if cds_feature.strand == "+":
            # Le RBS doit être entre [Start - 12] et [Start - 6]
            rbs_start_abs = max(0, cds_feature.start - upstream_dist_max) + 1  # 1-based
            rbs_end_abs = cds_feature.start - upstream_dist_min

        else:  # Brin -
            rbs_start_abs = cds_feature.start + upstream_dist_min
            rbs_end_abs = cds_feature.start + upstream_dist_max

        if rbs_start_abs > len(sequence) or rbs_end_abs < 1:
            return []

            # Ajustement des bornes pour la recherche
        search_start = max(0, rbs_start_abs - 1)
        search_end = min(len(sequence), rbs_end_abs)

        window_seq = sequence[search_start:search_end]

        for i in range(len(window_seq)):
            if window_seq[i : i + len(motif_pattern)] == motif_pattern:
                abs_pos = search_start + 1 + i  # 1-based

                rbs_list.append(
                    Feature(
                        type="RBS",
                        start=abs_pos,
                        end=abs_pos + len(motif_pattern) - 1,
                        strand="+",
                        sequence=motif_pattern,
                        length=len(motif_pattern),
                    )
                )

        return rbs_list

    @staticmethod
    def convert_coordinates(
        strand: str, start_1based: int, end_1based: int, total_length: int
    ) -> Tuple[int, int]:
        """
        Convertit les coordonnées 1-based trouvées sur un brin vers le génome global (+).
        Si strand est '+', retourne (start, end).
        Si strand est '-', retourne les coordonnées inversées sur le brin +.
        """
        if strand == "+":
            return start_1based, end_1based
        else:
            # Sur le brin -, l'index 0 correspond à la fin du brin +.
            # start/end sont en 1-based.
            s_0 = start_1based - 1
            e_0 = end_1based - 1

            # Conversion vers indices 0-based du brin +

            pos_start_plus_0 = (total_length - 1) - e_0
            pos_end_plus_0 = (total_length - 1) - s_0

            # Retour en 1-based
            return pos_start_plus_0 + 1, pos_end_plus_0 + 1


class AnnotationPanel(QWidget):
    """Panneau gauche contenant les modules de détection et d'annotation."""

    # Signaux pour communiquer avec le reste de l'application
    signal_cds_scan = pyqtSignal(int, bool, bool)  # (min_len, forward, reverse)
    signal_rbs_scan = pyqtSignal(str, str, str)  # (rbs_motif, promoter_35, promoter_10)
    signal_annotation_add = pyqtSignal(str, int, int)  # (feature_type, start, end)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumWidth(320)
        self.setup_ui()

    def setup_ui(self):
        main_layout = QVBoxLayout(self)

        # Création des onglets
        self.tabs = QTabWidget()

        # --- 1. Onglet CDS / ORF Finder ---
        tab_cds = QWidget()
        layout_cds = QFormLayout(tab_cds)

        grp_params = QGroupBox("Paramètres de détection")
        lay_grp = QFormLayout(grp_params)

        self.sp_min_len = QSpinBox()
        self.sp_min_len.setRange(50, 3000)
        self.sp_min_len.setValue(100)  # Par défaut > 99 nt
        self.sp_min_len.setSingleStep(10)

        self.cb_fwd = QCheckBox("Brin (+)")
        self.cb_fwd.setChecked(True)

        self.rev_cb = QCheckBox("Brin (-)")
        self.rev_cb.setChecked(True)

        cb_container = QWidget()
        cb_layout = QHBoxLayout(cb_container)
        cb_layout.setContentsMargins(0, 0, 0, 0)  # Supprime les marges internes
        cb_layout.addWidget(self.cb_fwd)
        cb_layout.addWidget(self.rev_cb)

        self.cmb_criteria = QComboBox()
        self.cmb_criteria.addItems(
            ["Cinit le plus éloigné de Cterm", "Le plus long par région"]
        )

        lay_grp.addRow("Longueur min (nt):", self.sp_min_len)
        lay_grp.addRow("Brins:", cb_container)
        lay_grp.addRow("Critère:", self.cmb_criteria)

        layout_cds.addRow(grp_params)

        btn_scan = QPushButton("Scanner la séquence")
        btn_scan.setStyleSheet(
            "background-color: #4CAF50; color: white; padding: 8px; font-weight: bold;"
        )
        btn_scan.clicked.connect(self.on_scan_clicked)
        layout_cds.addRow(btn_scan)

        self.tabs.addTab(tab_cds, "CDS / ORF")

        # --- 2. Onglet RBS & Promoteurs ---
        tab_rbs = QWidget()
        layout_rbs = QFormLayout(tab_rbs)

        grp_params = QGroupBox("Motifs et Paramètres")
        lay_grp = QFormLayout(grp_params)

        self.le_rbs_motif = QLineEdit("AGGAGGUAA")
        self.le_prom_35 = QLineEdit("TTGACA")
        self.le_prom_10 = QLineEdit("TATAAT")

        # Position relative (exemple simplifié pour l'UI)
        sp_range = QSpinBox()
        sp_range.setRange(5, 50)
        sp_range.setValue(35)

        lay_grp.addRow("Motif RBS:", self.le_rbs_motif)
        lay_grp.addRow("Boîte -35:", self.le_prom_35)
        lay_grp.addRow("Boîte -10:", self.le_prom_10)

        # TSS Management
        grp_tss = QGroupBox("Position du Start (TSS +1)")
        lay_tss = QHBoxLayout(grp_tss)
        self.btn_auto_tss = QPushButton("Auto-detect")
        self.lbl_manual_tss = QLabel("(ou cliquer sur la séquence)")
        lay_tss.addWidget(self.btn_auto_tss)
        lay_tss.addWidget(self.lbl_manual_tss, stretch=1)

        layout_rbs.addRow(grp_params)
        layout_rbs.addRow("Zone de recherche amont (nt):", sp_range)
        layout_rbs.addRow(grp_tss)

        btn_search = QPushButton("Rechercher motifs")
        btn_search.setStyleSheet(
            "background-color: #2196F3; color: white; padding: 8px; font-weight: bold;"
        )
        btn_search.clicked.connect(self.on_scan_clicked)
        layout_rbs.addRow(btn_search)

        self.tabs.addTab(tab_rbs, "RBS / Promoteurs")

        # --- 3. Onglet Annotation Manuelle ---
        tab_manual = QWidget()
        layout_manual = QVBoxLayout(tab_manual)

        grp_cursor = QGroupBox("Mode Curseur")
        lay_csr = QVBoxLayout(grp_cursor)
        self.lbl_cursor_info = QLabel("Position actuelle : - (nt)")
        self.lbl_cursor_info.setStyleSheet(
            "background-color: #333; color: lime; font-family: monospace; padding: 5px;"
        )

        btn_mode_toggle = QPushButton("Basculer Mode Édition")
        lay_csr.addWidget(self.lbl_cursor_info)
        lay_csr.addWidget(btn_mode_toggle, alignment=Qt.AlignmentFlag.AlignCenter)

        grp_add_feature = QGroupBox("Ajouter Feature")
        lay_feat = QFormLayout(grp_add_feature)

        self.cbx_type = QComboBox()
        self.cbx_type.addItems(["CDS", "Promoteur", "RBS", "Inconnu"])

        self.sp_feat_start = QSpinBox()
        self.sp_feat_end = QSpinBox()
        # On simule une plage, à connecter au chargement FASTA

        btn_add = QPushButton("➕ Ajouter à l'annotation")

        lay_feat.addRow("Type:", self.cbx_type)
        lay_feat.addRow("Start:", self.sp_feat_start)
        lay_feat.addRow("End:", self.sp_feat_end)
        lay_feat.addWidget(btn_add)

        layout_manual.addWidget(grp_cursor, stretch=1)
        layout_manual.addWidget(grp_add_feature)

        # Connexion slots (mock)
        btn_scan.clicked.connect(self.on_scan_clicked)
        btn_search.clicked.connect(self.on_scan_clicked)
        btn_add.clicked.connect(
            lambda: self.signal_annotation_add.emit("Manual", 0, 0)
        )  # À connecter aux inputs réels

        self.tabs.addTab(tab_manual, "Annotation")

        main_layout.addWidget(self.tabs)

    def update_tss_label(self, pos):
        """Mise à jour visuelle du Start Signal"""
        self.lbl_cursor_info.setText(f"Position actuelle : {pos} (nt)")

    def on_scan_clicked(self):
        # Émettre les signaux vers le backend
        sender = self.sender()

        if not isinstance(sender, QPushButton):
            return

        text = sender.text()

        if "Scanner la séquence" in text:
            # C'est le bouton CDS
            forward = self.cb_fwd.isChecked()
            reverse = self.rev_cb.isChecked()
            min_len = self.sp_min_len.value()
            print(f"[UI] Signal CDS émis: Min={min_len}, Fwd={forward}, Rev={reverse}")
            self.signal_cds_scan.emit(min_len, forward, reverse)

        elif "Rechercher motifs" in text:
            # C'est le bouton RBS/Promoteurs
            rbs_motif = self.le_rbs_motif.text()
            prom_35 = self.le_prom_35.text()
            prom_10 = self.le_prom_10.text()
            print(
                f"[UI] Signal RBS émis: Motif={rbs_motif}, P35={prom_35}, P10={prom_10}"
            )
            self.signal_rbs_scan.emit(rbs_motif, prom_35, prom_10)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Annotation Génomique")
        self.resize(1200, 700)

        # Widgets
        self.browser = SequenceBrowser()
        self.status_bar = QStatusBar()
        self.lbl_pos = QLabel("Prêt. Chargez un fichier FASTA.")

        # Panneau gauche
        self.left_panel = AnnotationPanel(self)

        # Connexions signaux
        self.browser.position_changed.connect(self.update_status_bar_cursor)

        # Connexion au backend
        self.left_panel.signal_cds_scan.connect(self.on_action_scan_cds)
        self.left_panel.signal_rbs_scan.connect(self.on_action_search_motifs)

        self.setup_ui()
        self.setup_menu()
        self.setStatusBar(self.status_bar)

    def setup_ui(self):
        central = QWidget()

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self.left_panel)
        splitter.addWidget(self.browser)
        splitter.setSizes([350, 850])

        main_layout = QVBoxLayout(central)
        main_layout.setContentsMargins(0, 0, 0, 0)

        main_layout.addWidget(splitter)
        self.setCentralWidget(central)

    def setup_menu(self):

        menu_bar = QMenuBar(self)
        self.setMenuBar(menu_bar)

        file_menu = menu_bar.addMenu("&Fichier")

        if not file_menu:
            return

        open_action = QAction("&Ouvrir FASTA", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self.open_fasta)
        file_menu.addAction(open_action)

    def open_fasta(self):
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Charger une séquence FASTA", "", "FASTA Files (*.fasta *.fa)"
        )
        if filepath:
            self.browser.load_fasta(filepath)

    def update_status_bar_cursor(self, pos, base):
        """Mise à jour du statusBar ET de l'input manuel dans le panneau gauche"""
        self.lbl_pos.setText(f"Position: {pos} nt | Base: {base}")
        self.status_bar.clearMessage()  # Nettoyer message précédent si présent

        # On met à jour le label du panneau gauche (Onglet 3) pour guider l'annotation
        self.left_panel.update_tss_label(pos)

    def on_action_scan_cds(self, min_len, fwd, rev):
        """Slot connecté au bouton 'Scanner' du panneau gauche"""
        print(f"--- LANCEMENT DÉTECTION CDS (Min={min_len}, Fwd={fwd}, Rev={rev}) ---")

        if not self.browser.sequence:
            self.status_bar.showMessage("Erreur : Aucune séquence chargée.", 5000)
            return

        seq_len = len(self.browser.sequence)
        cds_list = GenomeAnalyzer.find_orfs(
            sequence=self.browser.sequence, min_length=min_len, forward=fwd, reverse=rev
        )
        cds_list.sort(key=lambda x: x.start)
        print(f"  - Trouvé {len(cds_list)} CDS candidats.")
        msg = f"Détecté : {len(cds_list)} CDS. "
        for i, cds in enumerate(cds_list):
            start = cds.start
            end = cds.end

            if cds.strand == "-":
                new_start = seq_len - cds.end + 1
                new_end = seq_len - cds.start + 1
                print(
                    f"CDS #{i + 1}: Brin=-, Pos=[{new_start}..{new_end}], Len={cds.length}, StartATG={cds.sequence[:3]}"
                )
            else:
                print(
                    f"CDS #{i + 1}: Brin={cds.strand}, Pos=[{start}..{end}], Len={cds.length}"
                )
        self.status_bar.showMessage(msg, 5000)

    def on_action_search_motifs(self, rbs_motif, prom_35, prom_10):
        """Slot pour la recherche de motifs RBS et Promoteurs"""
        print(f"--- RECHERCHE MOTIFS ---")

        if not self.browser.sequence:
            return

        # Recherche RBS (ex: dans les 500 premiers nucléotides pour l'exemple)
        rbs_found = GenomeAnalyzer.find_motifs(
            sequence=self.browser.sequence,
            motif=rbs_motif,
            window_start=1,
            window_end=min(
                2000, len(self.browser.sequence)
            ),  # Limite arbitraire pour l'exemple
        )

        print(f"  - Motif RBS '{rbs_motif}' trouvé {len(rbs_found)} fois.")

    def update_status(self, pos, base):
        self.lbl_pos.setText(f"Position: {pos} nt | Base: {base}")


if __name__ == "__main__":
    app = QApplication(sys.argv)

    app.setStyle("Fusion")

    window = MainWindow()
    window.show()
    sys.exit(app.exec())
