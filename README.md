# circuit_simulation
Αλγόριθμοι Προσομοίωσης Κυκλωμάτων

Ακ.Έτος 2014 – 2015

Γλώσσα υλοποίησης : C

Το project περιλαμβάνει υλοποίηση ενός parser που διαβάζει netlist εισόδου με συγκεκριμένα χαρακτηριστικά όπως:
•	Non case sensitive
•	Μια γραμμή από κενά αντιλαμβάνονται ως ένα κενό
•	Κάθε γραμμή του netlist περιγράφει ένα μόνο κυκλωματικό στοιχείο
•	Γραμμή που αρχίζει με * περιέχει σχόλια και αγνοείται
•	Οι κόμβοι του κυκλώματος είναι ακέραιος ή αλφαριθμητικά (Χρήση hash Table)
•	O κόμβος 0(zero) αντιπροσωπεύει πάντα τον κόμβο αναφοράς ή γείωσης.
•	Οι τιμές των κυκλωματικών στοιχείων δίνονται χωρίς πολ/σια.
•	Κυκλωματικά στοιχεία που υποστηρίζονται είναι: 
Ανεξάρτητη πηγή τάσης:  V <name> <+> <-> <value>
Ανεξάρτητη πηγή ρεύματος: Ι <name> <+> <-> <value>
Αντίσταση: R <name> <+> <-> <value>
Χωρητικότητα: C <name> <+> <-> <value>
Αυτεπαγωγή: L <name> <+> <-> <value>
Δίοδος: D <name> <+> <-> <model-name> [area]
Τρανσίστορ MOS: M <name> <D> <G> <S> <B> <model-name>
Τρανσίστορ BJT: Q <name> <C> <B> <E>

Ο parser αποθηκεύει τα κυκλωματικά στοιχεία σε συνδεδεμένη λίστα όπου για κάθε στοιχείο-εγγραφή  καταγράφεται ο τύπος, οι ακροδέκτες και οι τιμή τους.
Ο κώδικας περιλαμβάνει χρήση Hash Table με τον οποίο διατρέχεται η συνδεδεμένη λίστα ( για οποιοδήποτε netlist εισόδου με στοιχεία RLC και πηγές) και δημιουργεί το DC σύστημα MNA. Για DC ανάλυση οι χωρητικότητες ανοιχοκυκλώνονται και οι αυτεπαγωγές βραχυκυκλώνονται (αντικαθίστανται από πηγές τάσης)
Υλοποιούνται Άμεσες(DIRECT) μέθοδοι επίλυσης Γραμμικών Συστημάτων, όπως Παραγοντοποίηση LU(Gauss) και Παραγοντοποίηση Cholesky χρησιμοποιούνται για την επίλυση συστήματος MNA. Η ένδειξη .OPTIONS SPD  στο netlist ενεργοποιεί την Cholesky ενώ η απουσία της ενεργοποιεί την LU. Εκτός της DC ανάλυσης προστέθηκε και η δυνατότητα σάρωσης DC (DC SWEEP). 
Επίσης υλοποιούνται Επαναληπτικές (Iterative) μέθοδοι επίλυσης γραμμικών συστημάτων όπως η Μέθοδος Συζυγών Κλίσεων ( Conjugate Gradients ) και η μέθοδος για μή συμμετρικά συστήματα με γρήγορο ρυθμό σύγκλισης Bi-CG.
Τέλος, υλοποιούνται οι παραπάνω μέθοδοι επίλυσης και για υποστήριξη Αραιών πινάκων(SPARSE). 

**Το Project έχει υλοποιηθεί με τρόπο που να υποστηρίζει netlist εισόδου μεγάλων κυκλωμάτων (IBM)
