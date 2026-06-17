# ☄️ Univerzális Meteorraj Tervező & Naptár Generáló

Egy intelligens Python script, amely csillagászati pontossággal számítja ki a négy legnagyobb meteorraj (**Lyridák, Perseidák, Orionidák, Geminidák**) láthatóságát Magyarország (alapértelmezetten Budapest) területére vetítve. 

A script nemcsak a meteorrajok elméleti csúcspontját figyeli, hanem **eleve a megfigyelési időablakok alatt vizsgálja a Hold pozícióját, horizont feletti magasságát és fázisát**. Ezzel pontosan megmutatja, hogy a fényszennyezés melyik éjszakákon fogja elnyomni a hullócsillagokat, és mikor kapunk tökéletesen sötét égboltot a kerti pokrócozáshoz.

---

## ✨ Főbb funkciók

* **4 nagy meteorraj támogatása:** Lyridák (április), Perseidák (augusztus), Orionidák (október) és Geminidák (december).
* **Egyedi időablakok:** A Geminidák korai keléséhez igazított 21:00-s indítás, míg a többi rajhoz a klasszikus 22:00-s sáv kezelése.
* **Intelligens Hold-analízis:** Kiszámítja a Hold pontos fázisát és horizont feletti magasságát, és ez alapján automatikusan text-alapú ajánlásokat/osztályzatokat generál.
* **Kettős kimenet:**
  1. 📄 **Konzol & `.txt` jelentés:** Egy szépen formázott áttekintő táblázat az elkövetkező évekről.
  2. 📅 **`.ics` naptárfájl:** Google Naptárba, Apple Calendarba vagy Home Assistant alá közvetlenül importálható naptáresemények (egésznapos emlékeztetőkkel és részletes leírást tartalmazó időzített blokkokkal).
* **Szabad státusz:** A naptárbejegyzések `TRANSPARENT` foglaltságúak, így nem blokkolják a mindennapi naptáradat, nem jelölnek téged "elfoglaltnak".

---

## 📸 Képernyőképek

### Szöveges áttekintő táblázat (Konzol & TXT kimenet)
![Konzol kimenet](https://github.com/bbkbarbar/shooting-stars-calendar/blob/main/result.png)

### Google Naptár integráció (A generált .ics fájl importálása)
![Google Naptár importálás](https://github.com/bbkbarbar/shooting-stars-calendar/raw/main/google_calendar_import.png)

---

## 🛠️ Követelmények és telepítés

A script futtatásához **Python 3.10+** és az `astropy` csillagászati könyvtár szükséges a hozzá tartozó `numpy` és `pytz` csomagokkal együtt.
   ```bash
   pip install astropy
   ```

1. Klónozd a repót:
   ```bash
   git clone [https://github.com/bbkbarbar/shooting-stars-calendar.git](https://github.com/bbkbarbar/shooting-stars-calendar.git)
   cd shooting-stars-calendar
2. Futtatás:
   ```bash
   python hullocsillag_idoszakok.py {kezdő_évszám} {generáltandó_évek_száma}
   ```
   2026 és 2027 eseményeinek generálásához:
   ```bash
   python hullocsillag_idoszakok.py 2026 2
   ```
