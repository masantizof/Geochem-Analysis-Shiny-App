# 📊 Multivariate Geochemical Analysis Dashboard

![R](https://img.shields.io/badge/Language-R-blue) ![Shiny](https://img.shields.io/badge/Framework-Shiny-purple) ![ML](https://img.shields.io/badge/ML-Random%20Forest-green)

## 📋 Project Overview
Developed for the **Colombian Geological Survey (SGC)**, this interactive web application democratizes access to advanced statistical analysis. It allows geologists to perform multivariate analysis on rock samples without needing to write code.

The tool integrates SQL databases with a reactive frontend to analyze spectral and geochemical signatures.

## 🧠 Advanced Capabilities
* **Machine Learning:** Implementation of **Random Forest** algorithms to classify geological samples based on their spectral properties.
* **Dimensionality Reduction:** Automated **Principal Component Analysis (PCA)** to identify patterns in high-dimensional chemical datasets.
* **Interactive Filtering:** Users can filter millions of data points by Department, District, or Lithology in real-time.

## 🛠️ Tech Stack
* **Language:** R.
* **Web Framework:** Shiny.
* **Stats/ML:** `FactoMineR`, `randomForest`, `corrplot`.
* **Database:** SQLite / SQL integration.
