nella prima riga del file input.dat ci può essere uno 0 o un 1 (da modificare MANUALMENTE). 0 significa che il codice parte dalla configurazione cristallina, mentre 1 significa che il codice riprende dalla configurazione finale della run precedente.

Parto direttamente dalla temperatura target, senza cominciare con temperature più alte.

Ricordarsi di eliminare i file con clean se faccio modifiche e voglio runnare il codice di nuovo.

QUELLO CHE CONVIENE FARE E' EQUILIBRARE IL SISTEMA CON 3k - 10k nstep e poi fare la misura che mi interessa con nstep = 30k, per avere 30 blocchi da 100 misure ciascuno