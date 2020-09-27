from glowingmeme.clients.clients import Clients
from protocols.protocol_7_2.reports import Program, Assembly

# load relevant clients
cipapi, cellbase, cva = Clients().get_all_clients()

## test CVA queries
cases_client = cva.cases()

cases_iterator = cases_client.get_cases(
    program=Program.rare_disease,
    assembly=Assembly.GRCh38,
    caseStatuses="ARCHIVED_POSITIVE",
    hasCanonicalTrio=True,
)
case_1 = next(cases_iterator)

case = cases_client.get_case(identifier="43868", version="1")

report_events_client = cva.report_events()
report_events_iterator = report_events_client.get_report_events(
    caseId="43868", caseVersion="1"
)
report_event_1 = next(report_events_iterator)

variants_cva_client = cva.variants()
# variant = variants_cva_client.get_variants(program=Program.rare_disease, assembly=Assembly.GRCh38,
#     caseStatuses="ARCHIVED_POSITIVE", hasCanonicalTrio=True, variantId='86d62eb935341eb357f9d3a4f1aaeae7')
# variant = variants_cva_client.get_variant_by_id('86d62eb935341eb357f9d3a4f1aaeae7')
variant = variants_cva_client.get_variant_by_id("53816dff6dd23a1e1c48d0d6767c51c2")

## test cipapi queries
interpretation_request = cipapi.get_case(case_id="43868", case_version=1)
