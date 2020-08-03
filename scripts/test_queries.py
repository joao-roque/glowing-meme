from glowingmeme.clients.clients import Clients
from protocols.protocol_7_2.reports import Program, Assembly

cipapi, cellbase, cva = Clients().get_all_clients()

cases_client = cva.cases()

cases_iterator = cases_client.get_cases(
    program=Program.rare_disease, assembly=Assembly.GRCh38,
    caseStatuses="ARCHIVED_POSITIVE", hasCanonicalTrio=True)
case_1 = next(cases_iterator)

report_events_client = cva.report_events()
report_events_iterator = report_events_client.get_report_events(
    caseId="10021", caseVersion="1")
report_event_1 = next(report_events_iterator)

variants_cva_client = cva.variants()
variant = variants_cva_client.get_variants(program=Program.rare_disease, assembly=Assembly.GRCh38,
    caseStatuses="ARCHIVED_POSITIVE", hasCanonicalTrio=True)
# variant = variants_cva_client.get_variant_by_id('df3f170a296a4bd03d7b5a43d9d35a14')

for var in variant:
    hello = variants_cva_client.get_variant_by_id(var.id)
    x = 1