# RDP Runner – Kullanım Kılavuzu (TR)

Bu uygulama, Relationship Diagramming Process (RDP) adımlarını görsel ve adım adım takip edilebilir şekilde çalıştırmanıza yardımcı olur. REL matrisi üzerinden ilişkileri tanımlar, TCR ve sıralama (Pi) oluşturur, ardından yerleşim (layout) fazında WPV (Weighted Placement Value) hesaplarına göre yerleştirme yapar.

## 1) Gereksinimler ve Kurulum
- Python 3.9+ (uygulama 3.11 ile test edildi)
- Tkinter (genelde Python ile birlikte gelir)

Kurulum gerektirmez. Proje klasöründe aşağıdaki komutu çalıştırmanız yeterlidir:

```bash
python RDP.py
```

## 2) Arayüz Genel Bakış
Uygulama 3 panelli bir yapıya sahiptir:

### Sol Panel (Kontroller)
- **Number of Departments (3–20):** Bölüm sayısını belirler, matris otomatik üretilir.
- **Randomize Relations:** REL matrisini rastgele doldurur (A/E/I/O/U/X).
- **Configure Values (A/E/I/O/U/X):** Seçili hücreleri belirli ilişki değerine ayarlamak için kısayol butonları.
- **Start Calculation:** Hesaplamayı başlatır, TCR ve seçim adımlarını hazırlar.
- **Next Step:** Seçim/yerleşim adımlarını tek tek ilerletir.
- **Jump to Layout:** Seçim fazını tamamlayıp yerleşim fazını hazırlar (yerleştirme yapmaz).
- **Run to End:** Tüm adımları otomatik tamamlar.
- **Reset Steps:** Matrisi korur, hesaplama/yerleşim ilerlemesini sıfırlar.
- **Reset All:** Tüm verileri temizler, matrisi U’ya döndürür.
- **Preview Report:** Raporu uygulama içinde önizler.
- **Export Report (.txt):** Raporu metin dosyasına dışa aktarır.
- **Guide / RDP Steps Info:** Kısa yardım pencereleri.
- **Log panel:** Uygulamanın işlem günlüğü.

### Orta Panel
- **REL Matrix:** İlişki matrisi (düzenlenebilir, kaydırılabilir).
- **Layout Visualization:** Yerleşim görselleştirmesi (kaydırılabilir).

### Sağ Panel
- **TCR Table:** Bölüm başına TCR değerleri.
- **Sequence (Pi):** Seçim fazı sonucu sıralama.
- **Placement Log & WPV:** Yerleşim adımları ve WPV kayıtları.

## 3) Temel İş Akışı (Önerilen)
1. **Bölüm sayısını** belirleyin.
2. REL matrisini **elle düzenleyin** veya **Randomize Relations** ile doldurun.
3. **Start Calculation** ile hesaplamayı başlatın.
4. **Next Step** ile adım adım ilerleyin veya **Run to End** ile tüm süreci tamamlayın.
5. Yerleşim fazında her adımda WPV adayları (gri “ghost” kutular) gösterilir. **Next Step** ile onaylayıp yerleştirmeyi tamamlayın.
6. **Preview Report** ile raporu görüntüleyin, **Export Report** ile dışa aktarın.

## 4) REL Değerleri ve Anlamları
Uygulama A/E/I/O/U/X ilişkilerini destekler. Varsayılan ağırlıklar uygulamada sabittir ve raporda yer alır. REL değerlerinin anlamı kurum/derse göre değişebilir; bu nedenle kendi standardınıza göre kullanın.

## 5) Seçim (Selection) Fazı
- **Start Calculation** ile TCR değerleri hesaplanır.
- Seçim fazı adımları otomatik hazırlanır.
- **Next Step** ile seçim adımlarını tek tek izleyebilirsiniz.
- **Jump to Layout**, seçim adımlarını tamamlar ve yerleşim fazına hazırlar (yerleştirme başlamaz).

## 6) Yerleşim (Layout) Fazı
- Yerleşim fazı **Selection** tamamlandıktan sonra başlar.
- Her adımda uygun aday noktalar ve **WPV** değerleri gösterilir.
- “Ghost” kutular, aday yerleşimlerin görsel önizlemesidir.
- **Next Step** ile en iyi aday seçilir ve departman yerleştirilir.

## 7) Raporlama
Rapor içeriği:
- Girdi parametreleri ve REL matrisi
- TCR tablosu
- Nihai Sequence (Pi)
- Seçim adımları
- Yerleşim adımları
- Nihai yerleşim koordinatları + ASCII grid

## 8) Sık Yapılan Hatalar / İpuçları
- **TCR veya Sequence görünmüyorsa:** önce **Start Calculation** tıklayın.
- **Matrisi düzenleyemiyorsanız:** seçim/yerleşim başlamış olabilir; **Reset Steps** ile matrisi kilitten çıkarın.
- **Log panelinde güncelleme yoksa:** Start Calculation sonrası Next Step veya Run to End kullanın.
- **Sütun genişlikleri ilk açılışta eşit değilse:** pencereyi yeniden boyutlandırarak sasha zorlamayı tetikleyin.

## 9) Sınırlamalar
- Departman sayısı 3–20 aralığındadır.
- Uygulama tek dosyalı Tkinter projesi olarak tasarlanmıştır; performans büyük matrislerde sınırlı olabilir.

## 10) Sorun Giderme (Hızlı)
- **Export Report** çalışmıyorsa: yazma izni olan bir klasör seçtiğinizden emin olun.
- **Layout boş görünüyorsa:** seçim adımları tamamlanmamış olabilir; Start Calculation + Next Step/Run to End kullanın.

---

## 11) SSS / FAQ

### 1) Start Calculation tıklamadan neden TCR/Sequence görünmüyor?
Uygulama hesaplamaları **Start Calculation** sonrası başlatır. Bu butona basılmadan TCR ve Sequence oluşturulmaz.

### 2) Matrisi neden düzenleyemiyorum?
Seçim veya yerleşim süreci başladıktan sonra matris kilitlenir. **Reset Steps** ile ilerlemeyi sıfırlayıp matrisi tekrar düzenleyebilirsiniz. Tam sıfırlama için **Reset All** kullanın.

### 3) Jump to Layout ne yapar?
Seçim fazını (TCR ve sıralama adımlarını) tamamlar ve yerleşim fazına geçmeye hazırlar. Ancak **departmanları yerleştirmez**. Yerleştirme için **Next Step** ile ilerlemek gerekir.

### 4) WPV nedir?
WPV (Weighted Placement Value), bir aday konuma yerleştirmenin ilişkisel ağırlığını ifade eder. Uygulama her aday konum için WPV hesaplar ve en uygun konumu seçer.

### 5) Ghost kutular neyi gösterir?
Ghost kutular, yerleşim adımındaki **aday konumların** görsel önizlemesidir. Üzerindeki değerler WPV skorlarıdır.

### 6) Run to End ile Next Step arasındaki fark nedir?
**Next Step** adım adım ilerletir ve ara sonuçları görmenizi sağlar. **Run to End** tüm adımları otomatik tamamlar.

### 7) Report Preview ve Export Report aynı mı?
İçerik aynıdır. **Preview Report** uygulama içi önizleme açar, **Export Report** metin dosyası olarak kaydeder.

### 8) Export Report dosyası nereye kaydedilir?
Kayıt konumu, açılan dosya diyalogunda seçtiğiniz klasördür. Yazma izni olan bir dizin seçmelisiniz.

### 9) “Logs not updating” uyarısı alıyorum, ne yapmalıyım?
Önce **Start Calculation**’a tıklayın. Ardından **Next Step** veya **Run to End** ile ilerleyin.

### 10) Departman sayısı sınırı var mı?
Evet, **3–20** aralığı desteklenir.

### 11) Seçim fazı tamamlanmadıysa yerleşim neden görünmüyor?
Yerleşim fazı için seçim fazı tamamlanmalıdır. **Start Calculation** sonrası **Next Step** ile ilerleyin veya **Run to End** kullanın.

### 12) Raporda hangi bilgiler yer alır?
REL matrisi, TCR, Sequence, seçim adımları, yerleşim adımları ve nihai layout bilgileri raporda bulunur.

### 13) Uygulama açılışta panel boyutları düzgün değilse?
Uygulama panel (pane) konumlarını yükleme sonrası tekrar ayarlar. Yine de sorun yaşarsanız pencereyi yeniden boyutlandırarak düzenlemeyi tetikleyin.
