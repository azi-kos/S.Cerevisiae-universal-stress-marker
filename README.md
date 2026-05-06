# S. cerevisiae — Universal Stress Markers

Seminar iz funkcijske genomike. Računalniška analiza podatkov o genski ekspresiji kvasovke *Saccharomyces cerevisiae* v treh različnih stresnih pogojih z namenom iskanja skupnih odzivnih genov.

---

## Kaj sploh počnemo?

Vsaka celica ima ~6000 genov. Vsak gen se lahko izraža bolj ali manj — koliko se izraža, merimo s **sekvenciranjem RNA (RNA-seq)**. Več RNA = gen je bolj aktiven.

V tem projektu smo vzeli kvasovko v 4 različnih okoljih:
- **Normalno (kontrola)** — standardni pogoji
- **Hipotonično** — premalo osmotskega tlaka v okolju (celica nabrekne)
- **Stradanje aminokislin** — brez gradnikov za proteine
- **Stradanje glukoze** — brez glavnega vira energije

Za vsak pogoj smo imeli ~10-50 **posameznih celic** (single-cell RNA-seq) in za vsako celico izmerili aktivnost vsakega od ~6000 genov.

**Cilj:** Kateri geni so spremenili aktivnost v *vseh treh* stresnih pogojih? Ti geni verjetno predstavljajo splošni stresni odgovor kvasovke — ne glede na tip stresa.

---

## Podatki

- **Vir:** GEO dataset [GSE201386](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201386)
- **Organizem:** *Saccharomyces cerevisiae* (pekovski kvas), sev S288C
- **Metoda:** Single-cell RNA-seq (ena celica = en vzorec)
- **Vzorcev skupaj:** 117 celic

| Pogoj | Serije | Število celic |
|---|---|---|
| Kontrola (izotonično) | A, E, G, Y700 | 38 |
| Hipotonično | C, F | 12 |
| Stradanje aminokislin | D | 20 |
| Stradanje glukoze | Y200, Y500 | 47 |

Surovi podatki so prav tako podani v repozitoriju.

---

## Metode

### 1. Sestava podatkovne matrike
Iz 117 ločenih datotek smo sestavili eno matriko: **~6000 genov × 117 celic**.

### 2. Diferencialna ekspresijska analiza (DESeq2)
DESeq2 je standardno orodje za primerjavo genske ekspresije med pogoji. Za vsak gen statistično testira: *"Je ta gen bolj ali manj aktiven v stresnih celicah v primerjavi s kontrolo?"*

Rezultat za vsak gen:
- **log2FC** — kako močno se je gen spremenil (log2FC = 1 pomeni 2× več, log2FC = -1 pomeni 2× manj)
- **padj** — statistična zanesljivost, korigirana za večkratno testiranje (padj < 0.05 = signifikantno)

Filter za signifikantnost: **padj < 0.05** in **|log2FC| > 1** (vsaj 2-kratna sprememba).

### 3. Iskanje skupnih genov
Signifikantne gene iz vsakega pogoja smo primerjali med seboj — iskali smo **presek** (gene ki so signifikantni v vseh treh pogojih hkrati).

### 4. GO enrichment analiza
GO (Gene Ontology) je standardiziran sistem kategorij bioloških funkcij. Preverili smo ali so skupni geni statistično prekomerno zastopani v določenih funkcionalnih kategorijah — npr. "mRNA procesiranje" ali "DNA popravljanje".

---

## Rezultati

| Primerjava | Signifikantni geni |
|---|---|
| Hipotonično vs. kontrola | 383 |
| Stradanje aminokislin vs. kontrola | 2487 |
| Stradanje glukoze vs. kontrola | 3056 |
| **Skupni presek (vsi 3 pogoji)** | **150** |

### Skupni stresni odgovor
150 genov je signifikantno spremenilo aktivnost v vseh treh stresnih pogojih. GO analiza nakazuje enrichment v:
- **mRNA procesiranje** — celica selektivno regulira katere RNA se procesirajo
- **DNA popravljanje in rekombinacija** — stres povzroča poškodbe DNA
- **Fermentacija in etanolni metabolizem** — preklop na alternativne vire energije
- **Spliceosom** — regulacija izrezovanja intronov iz RNA

### Ključni gen: ZNF1 (YFL052W)
Najmočneje reguliran gen v vseh treh pogojih (log2FC ≈ 30). ZNF1 je **transkripcijski faktor** ki regulira preklop iz fermentacije na respiracijo in odpornost na osmotski stres — verjetno osrednji regulator splošnega stresnega odgovora.

---
Transkripcijski faktorji (ZNF1, PUL4, HMS1) — celica aktivira regulatorne proteine ki nato kaskadno spremenijo ekspresijo stotih genov. To razloži zakaj vidimo tako veliko skupnih DEG-ov.

GAL10 — aktivacija alternativnih sladkornih poti ko primanjkuje glukoze ali aminokislin je smiselna.

MEC3 — DNA popravljalni checkpoint je aktiviran v vseh treh stresih, kar potrjuje GO rezultat o DNA popravljanju.

SME1 — jedrni del spliceosoma, potrjuje GO rezultat o spliceosomu in mRNA procesiranju.

HSP33 — posebej zanimiv ker je homolog humanega PARK7/DJ-1 ki je povezan s Parkinsonovo boleznijo. Tudi pri kvasovki je odgovoren za odpornost na oksidativni stres in dolgožvost celic.

SDS3 — histon deacetilaza, potrjuje GO rezultat o epigenetski regulaciji.

DUS3 — edini konsistentno down gen med top 10, modificira tRNA. Celica morda zmanjša tRNA modifikacije da upočasni translacijo med stresom.

---


## Kako ponoviti analizo


### Zagon
1. Prenesi podatke z GEO (GSE201386) in jih razpaki v `podatki/GSE201386_RAW/`
2. Odpri `analiza_deseq2.R` v RStudiu in poženi
3. Poženi `go_analiza.R`
4. Poženi `vizualizacije.R`

---

## Omejitve

- **Normalizacija:** DESeq2 median-of-ratios predpostavlja stabilno skupno količino RNA med vzorci. Pri stradanju glukoze to morda ne drži popolnoma — celice z manj energije proizvajajo manj RNA globalno.
- **GO enrichment:** 150 skupnih genov je premalo za zanesljivo statistično analizo po Benjamini-Hochberg korekciji. GO rezultati so eksploratorni.
- **Single-cell variabilnost:** Visoka biološka variabilnost med posameznimi celicami (posebej pri stradanju glukoze) zmanjšuje statistično moč.