import datetime
import sys
import warnings
import astropy.coordinates
from astropy.coordinates import EarthLocation, get_body, get_sun, NonRotationTransformationWarning
from astropy.utils.exceptions import AstropyWarning
from erfa.core import ErfaWarning
from astropy.time import Time
import astropy.units as u
import numpy as np
import pytz

warnings.filterwarnings('ignore', category=NonRotationTransformationWarning)
warnings.filterwarnings('ignore', category=AstropyWarning)
warnings.filterwarnings('ignore', category=ErfaWarning)

# ==================== GLOBÁLIS BEÁLLÍTÁSOK ====================
LOCATION_NAME = "Budapest"
LAT = 47.554179
LON = 19.001864

POSITION_NAME = "Hármas határhegy"
ELEVATION = 127     # balatons: 103-105
YEAR_COUNT = 7  # Hány évet vizsgáljon előre

# METEORRAJOZ KONFIGURÁCIÓJA
# Hónap, csúcsnap, megfigyelés kezdete, ablak hossza (óra), vizsgált napok környezete a csúcs körül
METEOR_SHOWERS = {
    "Lyridák": {
        "month": 4, "peak_day": 22, "start_hour": 22, "duration_hours": 5,
        "days_to_check": [20, 21, 22, 23], "emoji": "🌱"
    },
    "Perseidák": {
        "month": 8, "peak_day": 12, "start_hour": 22, "duration_hours": 5,
        "days_to_check": [10, 11, 12, 13], "emoji": "☄️"
    },
    "Orionidák": {
        "month": 10, "peak_day": 21, "start_hour": 22, "duration_hours": 5,
        "days_to_check": [19, 20, 21, 22], "emoji": "🏹"
    },
    "Geminidák": {
        "month": 12, "peak_day": 13, "start_hour": 21, "duration_hours": 6, # Korábban kezdődik és hosszabb!
        "days_to_check": [11, 12, 13, 14], "emoji": "♊"
    }
}

LINE_WIDTH = 155

# ===============================================================

class LogWriter:
    """Intelligens író, ami a konzolra és egy fájlba is szinkronban menti a kimenetet"""
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "w", encoding="utf-8")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()

    def close(self):
        self.log.close()

def get_moon_phase_and_alt(astropy_time, location):
    moon_gcrs = get_body('moon', astropy_time, location).transform_to(astropy.coordinates.GCRS(obstime=astropy_time))
    sun_gcrs = get_sun(astropy_time).transform_to(astropy.coordinates.GCRS(obstime=astropy_time))
    
    altaz_frame = astropy.coordinates.AltAz(obstime=astropy_time, location=location, pressure=1013 * u.hPa)
    moon_local = moon_gcrs.transform_to(altaz_frame)
    
    v_moon = np.array([moon_gcrs.cartesian.x.value, moon_gcrs.cartesian.y.value, moon_gcrs.cartesian.z.value])
    v_sun = np.array([sun_gcrs.cartesian.x.value, sun_gcrs.cartesian.y.value, sun_gcrs.cartesian.z.value])
    v_moon_to_sun = v_sun - v_moon
    
    dot_product = np.dot(-v_moon, v_moon_to_sun)
    norm_moon = np.linalg.norm(v_moon)
    norm_moon_to_sun = np.linalg.norm(v_moon_to_sun)
    
    cos_phase = np.clip(dot_product / (norm_moon * norm_moon_to_sun), -1.0, 1.0)
    return (1.0 + cos_phase) / 2.0 * 100.0, moon_local.alt.deg

def get_sunset_time(current_day, location, local_tz):
    base_time = datetime.datetime.combine(current_day, datetime.time(15, 0)) # December miatt korábbi bázis
    for m in range(0, 540):
        check_dt = base_time + datetime.timedelta(minutes=m)
        check_utc = local_tz.localize(check_dt).astimezone(pytz.utc)
        sun = get_sun(Time(check_utc))
        sun_altaz = sun.transform_to(astropy.coordinates.AltAz(obstime=Time(check_utc), location=location))
        if sun_altaz.alt.deg <= -0.833:
            return local_tz.localize(check_dt)
    return None

def find_moon_events_for_window(date_target, start_hour, duration_hours, location, local_tz):
    """Megkeresi az adott ablak környezetéhez tartozó holdkeltét vagy lementét"""
    start_local = datetime.datetime.combine(date_target, datetime.time(start_hour, 0)) - datetime.timedelta(hours=4)
    start_utc = local_tz.localize(start_local).astimezone(pytz.utc)
    
    rise_time = None
    set_time = None
    
    start_time = Time(start_utc)
    start_altaz_frame = astropy.coordinates.AltAz(obstime=start_time, location=location, pressure=1013 * u.hPa)
    start_moon = get_body('moon', start_time, location)
    last_alt = start_moon.transform_to(start_altaz_frame).alt.deg + 0.25
    
    # Pásztázási ablak hossza: ablak hossza + 8 óra puffer (percekben)
    total_minutes = (duration_hours + 8) * 60
    
    for minute in range(5, total_minutes, 5):
        check_time_utc = start_utc + datetime.timedelta(minutes=minute)
        astropy_time = Time(check_time_utc)
        
        altaz_frame = astropy.coordinates.AltAz(obstime=astropy_time, location=location, pressure=1013 * u.hPa)
        moon = get_body('moon', astropy_time, location)
        current_alt = moon.transform_to(altaz_frame).alt.deg + 0.25
        
        local_dt = check_time_utc.astimezone(local_tz)
        
        if last_alt < 0 and current_alt >= 0 and not rise_time:
            rise_time = local_dt
        if last_alt >= 0 and current_alt < 0 and not set_time:
            set_time = local_dt
                
        last_alt = current_alt
        
    return rise_time, set_time

def check_shower_night(date_target, start_hour, duration_hours, location, local_tz):
    """Kiszámítja a Hold fentlétét a specifikus ablak alatt"""
    start_local = datetime.datetime.combine(date_target, datetime.time(start_hour, 0))
    moon_present_minutes = 0
    max_illumination = 0.0
    
    total_steps = (duration_hours * 60) + 10
    for m in range(0, total_steps, 10):
        check_dt = start_local + datetime.timedelta(minutes=m)
        check_utc = local_tz.localize(check_dt).astimezone(pytz.utc)
        illumination, alt = get_moon_phase_and_alt(Time(check_utc), location)
        if alt + 0.25 >= 0:
            moon_present_minutes += 10
            if illumination > max_illumination:
                max_illumination = illumination
                
    return max_illumination, moon_present_minutes / 60.0

if __name__ == "__main__":
    location = EarthLocation(lat=LAT * u.deg, lon=LON * u.deg, height=ELEVATION * u.m)
    local_tz = pytz.timezone("Europe/Budapest")
    
    current_year = datetime.date.today().year
    if len(sys.argv) > 1:
        try: current_year = int(sys.argv[1])
        except ValueError: pass

    if len(sys.argv) > 2:
        try: YEAR_COUNT = int(sys.argv[2])
        except ValueError: pass

    # Fájlnevek előkészítése
    txt_filename = f"hullocsillag_attekinto_{LOCATION_NAME.lower()}.txt"
    ics_filename = f"hullocsillag_programok_{LOCATION_NAME.lower()}.ics"

    # Átirányítjuk a rendszert, hogy a konzol mellett a TXT fájlba is írjon egyszerre
    logger = LogWriter(txt_filename)
    sys.stdout = logger

    ics_lines = [
        "BEGIN:VCALENDAR",
        "VERSION:2.0",
        f"PRODID:-//{LOCATION_NAME} Csillagaszati Naptar//NONSGML HU",
        "CALSCALE:GREGORIAN"
    ]

    print(f"✨ METEORRAJ-TERVEZŐ ({LOCATION_NAME}) ✨")
    print(f"Vizsgált időszak: {current_year} ---> {current_year + YEAR_COUNT - 1}")

    for yr in range(current_year, current_year + YEAR_COUNT):
        for shower_name, info in METEOR_SHOWERS.items():
            emo = info["emoji"]
            sh_month = info["month"]
            sh_peak = info["peak_day"]
            st_hour = info["start_hour"]
            dur_hours = info["duration_hours"]
            
            print("\n" + "="*LINE_WIDTH)
            print(f"{emo} {shower_name.upper()} METEORRAJ - {yr} OSZTÁLYOZÁS")
            print("-"*LINE_WIDTH)
            print(f"{'Dátum / Éjszaka':<15} | {'Napnyugta':<9} | {'Holdkelte':<9} | {'Holdlemente':<11} | {'Hold fázis':<10} | {'Hold fentlét':<12} | Ajánlás")
            print("-"*LINE_WIDTH)
            
            # 1. EGÉSNAPOS ESEMÉNY A CSÚCSRÓL AZ ICS-BE
            all_day_start = f"{yr}{sh_month:02d}{sh_peak:02d}"
            # Következő nap kiszámítása az egésznapos lezáráshoz
            peak_date = datetime.date(yr, sh_month, sh_peak)
            next_day_date = peak_date + datetime.timedelta(days=1)
            all_day_end = next_day_date.strftime("%Y%m%d")
            
            sunset_max = get_sunset_time(peak_date, location, local_tz)
            m_rise_max, m_set_max = find_moon_events_for_window(peak_date, st_hour, dur_hours, location, local_tz)
            
            m_rise_max_str = m_rise_max.strftime('%H:%M') if m_rise_max else "nem kel"
            m_set_max_str = m_set_max.strftime('%H:%M') if m_set_max else "nem nyugszik"
            sunset_max_str = sunset_max.strftime('%H:%M') if sunset_max else "--:--"

            all_day_desc = (
                f"A(z) {shower_name} meteorraj elméleti csúcspontja. A legintenzívebb hullócsillag-zápor az esti megfigyelési ablakban várható.\\n\\n"
                f"Csillagászati menetrend mára:\\n"
                f"• Napnyugta: {sunset_max_str}\\n"
                f"• Holdkelte: {m_rise_max_str}\\n"
                f"• Holdlemente: {m_set_max_str}"
            )

            ics_lines.extend([
                "BEGIN:VEVENT",
                f"DTSTART;VALUE=DATE:{all_day_start}",
                f"DTEND;VALUE=DATE:{all_day_end}",
                f"SUMMARY:{emo} {shower_name} - MAXIMUMA",
                f"DESCRIPTION:{all_day_desc}",
                fr"LOCATION:{LOCATION_NAME}\, Magyarország",
                "TRANSP:TRANSPARENT",
                "END:VEVENT"
            ])

            # 2. RÉSZLETES NAPI BONTÁS (Konzol + TXT + Időzített ICS)
            for day_num in info["days_to_check"]:
                target_date = datetime.date(yr, sh_month, day_num)
                max_ill, moon_hrs = check_shower_night(target_date, st_hour, dur_hours, location, local_tz)
                sunset_dt = get_sunset_time(target_date, location, local_tz)
                moon_rise, moon_set = find_moon_events_for_window(target_date, st_hour, dur_hours, location, local_tz)
                
                day_str = target_date.strftime('%Y-%m-%d')
                sunset_str = sunset_dt.strftime('%H:%M') if sunset_dt else "--:--"
                m_rise_str = moon_rise.strftime('%H:%M') if moon_rise else "nem kel"
                m_set_str = moon_set.strftime('%H:%M') if moon_set else "nem nyugszik"
                    
                if moon_hrs == 0: 
                    rec = "💎 TÖKÉLETES! Sötét égbolt végig, a Hold nem lesz fent az ablakban!"
                    max_ill = 0.0
                elif max_ill < 25: 
                    rec = "👍 KIVÁLÓ! Csak egy vékony sarló van fent, alig zavar."
                elif moon_hrs <= 2.0: 
                    rec = "⏳ JÓ! A Holdnak csak rövid szakasza van fent, van tiszta, sötét ablak."
                elif max_ill > 80: 
                    rec = "🛑 ROSSZ! Erős holdfény, elnyomja a hullócsillagokat."
                else: 
                    rec = "⚠️ KÖZEPES! A Hold zavarni fogja az észlelést."
                    
                # Ez a sor a LogWriter-nek köszönhetően a képernyőre és a TXT-be is beíródik egyszerre!
                print(f"{day_str:<15} | {sunset_str:<9} | {m_rise_str:<9} | {m_set_str:<11} | {max_ill:>9.1f}% | {moon_hrs:>8.1f} óra | {rec}")

                # Időzített egyedi bejegyzések az ICS-be
                start_local = datetime.datetime.combine(target_date, datetime.time(st_hour, 0))
                end_local = start_local + datetime.timedelta(hours=dur_hours)
                
                start_utc = local_tz.localize(start_local).astimezone(pytz.utc).strftime("%Y%m%dT%H%M%SZ")
                end_utc = local_tz.localize(end_local).astimezone(pytz.utc).strftime("%Y%m%dT%H%M%SZ")
                
                summary = f"{emo} Hullócsillag nézés ({shower_name})"
                if day_num == sh_peak: summary += " [CSÚCS éjszaka]"
                
                desc = (
                    f"Megfigyelési ablak: {start_local.strftime('%H:%M')} - {end_local.strftime('%H:%M')}\\n\\n"
                    f"Pontos adatok az éjszaka környezetében:\\n"
                    f"• Napnyugta: {sunset_str}\\n"
                    f"• Holdkelte: {m_rise_str}\\n"
                    f"• Holdlemente: {m_set_str}\\n"
                    f"• Hold max telítettsége az ablak alatt: {max_ill:.1f}%\\n"
                    f"• Hold fentlét a megfigyelési időben: {moon_hrs:.1f} óra\\n\\n"
                    f"Ajánlás: {rec}"
                )

                ics_lines.extend([
                    "BEGIN:VEVENT",
                    f"DTSTART:{start_utc}",
                    f"DTEND:{end_utc}",
                    f"SUMMARY:{summary}",
                    f"DESCRIPTION:{desc}",
                    fr"LOCATION:{LOCATION_NAME}\, {POSITION_NAME}",
                    "TRANSP:TRANSPARENT",
                    "END:VEVENT"
                ])
            print("="*LINE_WIDTH)

    ics_lines.append("END:VCALENDAR")
    
    with open(ics_filename, "w", encoding="utf-8") as f:
        f.write("\n".join(ics_lines))
        
    # Visszaállítjuk az eredeti kimenetet, hogy a sikeres mentés üzenet már csak a konzolon fusson le
    sys.stdout = logger.terminal
    logger.close()
        
    print(f"\n✨ Sikeres mentés!")
    print(f" 📂 Összesített naptárfájl létrehozva: {ics_filename}")
    print(f" 📄 Szöveges áttekintő táblázat elmentve: {txt_filename}")