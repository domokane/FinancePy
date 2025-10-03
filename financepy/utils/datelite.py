from .date import Date

class DateLite:
    __slots__ = ("excel_dt",)

    def __init__(self, excel_dt: int):
        self.excel_dt = int(excel_dt)

    # Add days quickly
    def add_days(self, n: int) -> "DateLite":
        return DateLite(self.excel_dt + n)

    # Weekday using Excel convention (0=Mon … 6=Sun)
    def weekday(self) -> int:
        return (self.excel_dt + 5) % 7

    def is_weekend(self) -> bool:
        wd = self.weekday()
        return wd == 5 or wd == 6

    # Convert back to full Date when needed
    def to_date(self) -> "Date":
        return Date.from_excel(self.excel_dt)

    def __repr__(self):
        return f"DateLite({self.excel_dt})"

