using IsoconversionalAnalysis
using Logging

project_root = normpath(joinpath(@__DIR__, ".."))
catalog_path =
    isempty(ARGS) ? joinpath(project_root, "config", "datasets.toml") : abspath(ARGS[1])
markdown_path =
    length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(project_root, "docs", "data_audit.md")
toml_path = if length(ARGS) >= 3
    abspath(ARGS[3])
else
    joinpath(project_root, "docs", "audits", "m2_data_inventory.toml")
end

catalog = load_dataset_catalog(catalog_path)
@info "data_audit_started" catalog = catalog.source_path
audits = audit_catalog(catalog)
outputs = write_audit_reports(
    catalog, audits; markdown_path=markdown_path, toml_path=toml_path
)
@info "data_audit_completed" experiment_count = length(audits) markdown = outputs.markdown toml =
    outputs.toml
