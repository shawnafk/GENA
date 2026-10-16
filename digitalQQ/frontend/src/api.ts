export async function request<T>(
  url: string,
  data?: unknown,
  signal?: AbortSignal,
): Promise<T> {
  const r = await fetch(url, {
    method: data === undefined ? "GET" : "POST",
    headers: { "Content-Type": "application/json" },
    body: data === undefined ? undefined : JSON.stringify(data),
    signal,
  });
  if (!r.ok) {
    const e = await r.json();
    throw new Error(
      typeof e.detail === "string" ? e.detail : JSON.stringify(e.detail),
    );
  }
  return r.json();
}
// A generation gate rejects stale successes AND stale errors even when fetch abort is too late.
export class Latest {
  private generation = 0;
  next() {
    return ++this.generation;
  }
  current(id: number) {
    return id === this.generation;
  }
}
export function save(
  name: string,
  content: BlobPart,
  type = "application/json",
) {
  const url = URL.createObjectURL(new Blob([content], { type }));
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
export function localTime(iso: string) {
  return new Date(new Date(iso).getTime() + 8 * 3600000)
    .toISOString()
    .slice(0, 16);
}
export function utc(local: string) {
  return new Date(local + "+08:00").toISOString();
}
