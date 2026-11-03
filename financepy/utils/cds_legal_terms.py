##############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
##############################################################################


from dataclasses import dataclass


@dataclass
class CDSLegalTerms:
    reference_entity: str = ""
    reference_obligation: str = ""
    currency: str = ""
    seniority: str = ""
    restructuring: str = ""
    documentation: str = ""
    settlement_method: str = ""
