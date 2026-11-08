# Сравнительный анализ поверхностно-аннотированных белков Candida auris

Snakemake-пайплайн проекта ТЕН. Он объединяет каталог UniProtKB, кластеры сходства MMseqs2, поиск DIAMOND, контроль качества сборок, ANI/AF, картирование miniprot, выравнивание MAFFT и оценку консервативности. Итоговый отчёт сопоставляет вычислительные свидетельства и ограничения для каждого белка.

## Быстрый старт

Проверенная платформа запуска — Linux x86-64, CPU. Команды выполняются из корня репозитория. Для реальных данных укажите пути в копии конфигурации.

```bash
conda env create -f environment.yml
conda activate cauris-workflow
conda config --env --set channel_priority strict
cp config/config.yaml config/local.yaml
# Укажите пути к своим входным файлам в config/local.yaml
snakemake --cores 16 --configfile config/local.yaml --dry-run
snakemake --cores 16 --sdm conda --rerun-incomplete --show-failed-logs \
  --configfile config/local.yaml
```

Snakemake создаёт отдельные окружения для групп правил. Версии основных инструментов зафиксированы в [environment.yml](environment.yml) и [workflow/envs/](workflow/envs/). Для повторяемого запуска используйте одинаковые входные файлы, конфигурацию и окружения. JSON-метаданные и ZIP могут различаться между прогонами из-за времени и путей.

### Проверка установки без исходного набора

```bash
python tests/make_fixture.py
snakemake all legacy_mapping_summary --cores 4 --sdm conda \
  --show-failed-logs --configfile tests/_work/config.json
python tests/check_fixture.py tests/_work/results --legacy
```

Тест генерирует искусственные последовательности. Он проверяет чистое извлечение белков, несколько аллелей, единственный аллель, отсутствие выравнивания, DIAMOND, skani, SHA-256 и содержимое итогового ZIP. [GitHub Actions](.github/workflows/smoke.yml) запускает эту проверку, затем проверяет повторный запуск и режим архивных матриц ANI/AF.

## Входные данные

| Ключ в конфиге | Что подготовить |
| --- | --- |
| `catalog.directory` | Снимок UniProtKB с файлами `surface_candidates.faa`, `surface_candidates.tsv` и `run_metadata.json`. При наличии также копируются `query_audit.tsv`, `by_organism.tsv` и `cell_wall_annotated.faa`. |
| `inputs.selected` | Список accession сборок с версиями, по одному на строку. Поддерживается исходный TSV со столбцом `Assembly Accession`. Для расчёта ANI требуется минимум две сборки. |
| `inputs.genomes` | Каталог NCBI Datasets `ncbi_dataset/data/`. В папке каждого выбранного accession должен быть ровно один геномный FASTA, обычный или gzip. |
| `inputs.human_proteome` | FASTA референсного протеома человека. Для исторического запуска возьмите тот же выпуск и тот же набор изоформ. |
| `inputs.clade_assignments` | Архивная таблица с колонками `assembly`, `assigned_clade`, `assignment_status` и опциональной `known_control_clade`. |
| `inputs.clade_comparisons` | Таблица сравнений с референсами с колонками `assembly`, `ref_clade`, `reference_accession`, `ANI_pct`, `AF_ref_pct`, `AF_query_pct`. |
| `inputs.candidate_accessions` | Необязательный список accession белков. Позволяет зафиксировать исторический набор при наличии исходного снимка каталога. |
| `inputs.isolates` | Необязательный каталог с файлами `<assembly>/protein.faa` для отдельного анализа аннотированных изолятов. |
| `inputs.ani_matrix`, `inputs.af_matrix` | Необязательная пара сохранённых матриц. Если обе отсутствуют, skani пересчитает матрицы. AF должен содержать оба направления. |

Для текстового списка accession геномы можно получить через NCBI Datasets. Сохраняйте версии сборок и исходный архив вместе с метаданными запуска. [Документация NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/)

```bash
datasets download genome accession --inputfile inputs/selected_accessions.txt \
  --include genome --filename inputs/ncbi_genomes.zip
unzip inputs/ncbi_genomes.zip -d inputs
```

Для нового каталога можно переключить `catalog.mode` на `live`. Пайплайн выполнит запросы к UniProtKB, сохранит аудит запросов и идентификатор выпуска, если API его возвращает. Чтобы повторить полученный анализ, перенесите экспорт из `results/01_surface_antigens/data/` в отдельный каталог входных данных и используйте режим `snapshot`. [Документация UniProt](https://www.uniprot.org/help/api)

## Результаты

| Каталог под `results/` | Содержание |
| --- | --- |
| `00_provenance/` | SHA-256 входов и файлов workflow, итоговая конфигурация |
| `02_qc_catalog/`, `03_clusters/`, `04_summarize_clusters/` | QC каталога, группы сходства, FASTA приоритетной панели |
| `06_blastp_on_human/` | Хиты DIAMOND и сводка для всех кандидатов, включая записи без хитов |
| `07_annotated_isolates/summary/` | Отдельная сводка по аннотированным протеомам при включённой ветке |
| `08_population/` | QC сборок, ANI/AF, диагностика клад, картирование и аудит |
| `09_protein_sequences/`, `10_mafft/`, `11_conservation_statistics/` | Извлечённые последовательности, аллели, выравнивания, статистика по позициям и кладам |
| `12_final_report/` | Итоговые TSV, русский отчёт, графики, манифесты и `final_results_bundle.zip` |
| `13_legacy_population/` | Историческая сводка по правилам первой версии популяционного скрипта |

Для сравнения исторической и основной сводки добавьте цель `legacy_mapping_summary` к запуску, как в примере теста. Логи правил хранятся в `results/logs/`, а логи отдельных выравниваний также в каталогах соответствующих этапов.

## Воспроизведение исходного анализа и публикация

В приложенных скриптах не было кода первоначального отбора сборок и присвоения клад. Эти решения оформлены как версионируемые входные данные. Таблица ANI всей панели используется для аудита сходства, а метки клад берутся из отдельной таблицы назначений. В [заметках по переносу](docs/migration.md) указаны правила, добавленные между исходными скриптами, и параметры, которые следует сверить с журналами исходного прогона.

Перед публичным релизом укажите авторов, выберите лицензию на код и добавьте ссылки на архивные снимки входов с их контрольными суммами. Исходные данные, локальная конфигурация и результаты исключены из индекса Git через `.gitignore`.

Группы MMseqs2 описывают сходство последовательностей. Расчёты консервативности относятся к успешно извлечённым гомологам. Для выводов об экспрессии, поверхностной локализации, перекрёстной реактивности и иммунной защите нужны независимые данные.
