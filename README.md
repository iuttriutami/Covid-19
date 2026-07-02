# COVID-19 Cluster Analysis — Shiny App

Aplikasi **R Shiny** interaktif untuk analisis klaster COVID-19 34 provinsi
Indonesia, dikonversi dari script `covid19_cluster_analysis_COMPLETE_EN_v3.R`
(model manuskrip IJDNS).

🔗 Data: kasus/100k penduduk, CFR, tingkat kesembuhan, cakupan vaksinasi dosis-2
(kumulatif s.d. Des 2022). Data sudah ter-embed di `app.R`, jadi aplikasi bisa
langsung jalan tanpa file eksternal.

## Fitur

| Tab | Isi |
|---|---|
| **Data** | Tabel data mentah 34 provinsi + download CSV |
| **Preprocessing & K** | Pilih transformasi (log), skala (z-score/min-max/robust), dosis vaksin; plot Elbow/Silhouette/Gap statistic |
| **Clustering & Validation** | Jalankan 4 metode (K-Means, Ward, GMM, Fuzzy C-Means), tabel validasi (Silhouette, Davies-Bouldin, Calinski-Harabasz, Dunn), metode terbaik otomatis (voting), profil klaster |
| **Bootstrap Stability** | Stabilitas klaster via Jaccard bootstrap (jumlah resample bisa diatur) |
| **Sensitivity Analysis** | A) jumlah k, B) preprocessing, C) dosis vaksin, D) metrik jarak Ward, E) konkordansi antar-metode (ARI), F) fuzzifier FCM, G) leave-one-province-out |
| **Visualizations** | PCA 4 metode, Silhouette plot, Dendrogram Ward |
| **Cluster Map** | Peta choropleth — upload shapefile (.shp/.shx/.dbf/.prj) sendiri, join via kode provinsi BPS (kolom `PROVNO`) |
| **Export** | Unduh hasil klaster, validasi, dan stabilitas bootstrap sebagai CSV |

## Menjalankan secara lokal

```r
# 1. Install dependencies
install.packages(c(
  "shiny", "shinythemes", "DT", "ggplot2", "dplyr", "tidyr",
  "cluster", "clusterSim", "fpc", "factoextra", "mclust", "e1071",
  "dendextend", "ggrepel", "patchwork", "knitr", "RColorBrewer",
  "sf", "shinycssloaders"
))

# 2. Jalankan aplikasi
shiny::runApp("app.R")
```

Paket `clusterSim`, `fpc`, dan `sf` bersifat opsional — jika tidak
terinstal, aplikasi tetap berjalan tapi beberapa indeks validasi (Davies-Bouldin,
Calinski-Harabasz, Dunn) atau tab peta akan dinonaktifkan/menampilkan `NA`.

## Deploy ke shinyapps.io

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name = "<akun-anda>",
                           token = "<token>",
                           secret = "<secret>")
rsconnect::deployApp(appDir = ".")
```

## Deploy via Docker (opsional)

```dockerfile
FROM rocker/shiny:4.4.0
RUN R -e "install.packages(c('shinythemes','DT','ggplot2','dplyr','tidyr', \
  'cluster','clusterSim','fpc','factoextra','mclust','e1071','dendextend', \
  'ggrepel','patchwork','knitr','RColorBrewer','sf','shinycssloaders'), \
  repos='https://cran.rstudio.com/')"
COPY app.R /srv/shiny-server/covid19-cluster/app.R
EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
```

## Upload ke GitHub

```bash
git init
git add app.R README.md .gitignore
git commit -m "Add COVID-19 cluster analysis Shiny app"
git branch -M main
git remote add origin https://github.com/<username>/<repo-name>.git
git push -u origin main
```

Setelah itu, repo bisa langsung dihubungkan ke **shinyapps.io** atau
**Posit Connect** untuk hosting otomatis, atau dijalankan siapa saja dengan
`shiny::runGitHub("<username>/<repo-name>")`.

## Struktur file

```
.
├── app.R          # Aplikasi Shiny (UI + server, satu file)
├── README.md      # Dokumen ini
└── .gitignore
```

## Catatan peta (Cluster Map)

Shapefile batas provinsi (`indo_by_prov_2023.shp` dkk.) **tidak** disertakan
di repo ini agar ukurannya tetap kecil dan portable. Upload file shapefile
Anda sendiri di tab **Cluster Map** (harus punya kolom `PROVNO` berisi kode
wilayah BPS 2 digit, sama seperti kolom `code` di data provinsi).

## Referensi metode

- Kaufman & Rousseeuw (1990) — Silhouette
- Davies & Bouldin (1979) — Davies-Bouldin Index
- Calinski & Harabasz (1974) — Calinski-Harabasz Index
- Hubert & Arabie (1985) — Adjusted Rand Index
- Fang & Wang (2012) — Bootstrap cluster stability
