Vorwort:
Das beiliegende Programm entstand als Teil meiner Bachelorarbeit.
Es existiert eine alte Version, die eine Inkonsistenz im sweep des Monte-Carlo-Verfahrens hat,
womit ich die Ergebnisse aus der Bachelorarbeit erzeugt habe.
Hierbei wird zwischen jedem Wert ein sweep mehr durchgeführt als geplant. Dies hat jedoch keinen Einfluss
auf die Aussagekraft der Ergebnisse. In der Standardversion ist dieser Fehler behoben.

Angaben zum Programm:
Das Programm verwendet den Metropolis sweep und das Monte-Carlo-Verfahren für die Simulation eines SU(2)-Gitters.
Hierbei wird der Fehler über bootstrap copies berechnet und sowohl eine Modifizierung des Metropolis Sweep über
das staple, als auch die binning Methode verwendet.
Es handelt sich um einen linearen Algorithmus, der Schritt für Schritt alles abarbeitet. Somit ist das Programm
nicht für komplexe Simulationen oder großen Gittervolumina geeignet.


Copyright:TBA

Kompilierung:
g++ -std=c++14 Monte_Carlo_SU2.cpp -o Monte_Carlo_SU2
./Monte_Carlo_SU2