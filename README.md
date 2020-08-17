Autor: Tobias Grundmann

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


Copyright: Creative Commons Lizenz vom Typ Namensnennung 3.0 Deutschland

Kompilierung:
g++ -std=c++14 Monte_Carlo_SU2.cpp -o Monte_Carlo_SU2
./Monte_Carlo_SU2

Im Ordner Bilder befinden sich die verwendeten Bilder für die Bachelorarbeit.
Im Ordner Daten befinden sich die Ausgaben des Programms.
Im Ordner Ergebnisse befinden sich die Werte der Datenpunkte aus den Bildern,
die Namen der Datei und des Bildes sind bis auf die Dateiendung identisch.