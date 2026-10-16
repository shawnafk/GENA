import { chromium } from "@playwright/test";
import { mkdirSync } from "node:fs";
import { resolve } from "node:path";
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1500, height: 1050 } });
const errors = [];
page.on("pageerror", (e) => errors.push(e.message));
await page.goto("http://127.0.0.1:8765");
await page
  .getByTestId("scene-plot")
  .locator("canvas")
  .first()
  .waitFor({ timeout: 60000 });
await page.getByText("目标可达", { exact: true }).waitFor({ timeout: 120000 });
await page.waitForFunction(
  () =>
    document.querySelector(".header-meta")?.textContent?.includes("场景就绪"),
  {},
  { timeout: 120000 },
);
mkdirSync(resolve("../artifacts"), { recursive: true });
await page.screenshot({
  path: resolve("../artifacts/desktop.png"),
  fullPage: true,
});
const imageDownload = page.waitForEvent("download");
await page.getByRole("button", { name: "导出当前 PNG", exact: true }).click();
await (await imageDownload).saveAs(resolve("../artifacts/browser_scene.png"));
await page.getByLabel("回放起点").fill("2026-10-21T00:00");
await page.getByLabel("回放终点").fill("2026-10-21T01:00");
const reportDownload = page.waitForEvent("download", { timeout: 120000 });
await page
  .getByRole("button", { name: "导出观测报告 ZIP", exact: true })
  .click();
await (await reportDownload).saveAs(resolve("../artifacts/browser_report.zip"));
await page.getByLabel("回放起点").fill("2026-10-20T00:00");
await page.getByLabel("回放终点").fill("2026-10-22T00:00");
await page.getByLabel("theta_roll · °").fill("-25");
await page.waitForFunction(
  () =>
    document.querySelector(".header-meta")?.textContent?.includes("场景就绪"),
  {},
  { timeout: 60000 },
);
await page.getByRole("button", { name: "反解 roll / pitch →" }).click();
await page.waitForFunction(
  () =>
    Math.abs(
      Number(
        document.querySelector('input[aria-label="theta_roll · °"]').value,
      ) + 30.932074,
    ) < 1e-5,
);
await page.getByRole("button", { name: "XY", exact: true }).click();
await page.getByLabel("目标定义").selectOption("far");
await page.getByRole("button", { name: "反解 roll / pitch →" }).click();
await page.getByText("存在残差", { exact: true }).waitFor({ timeout: 60000 });
await page.getByLabel("姿态保持方式").selectOption("inertial");
await page
  .getByLabel("时间轴")
  .fill(String(new Date("2026-10-21T04:00:00Z").getTime()));
await page.waitForFunction(
  () =>
    document.querySelector(".header-meta")?.textContent?.includes("场景就绪"),
  {},
  { timeout: 120000 },
);
await page.getByLabel("姿态保持方式").selectOption("tracking");
await page.getByLabel("目标定义").selectOption("point");
await page.getByText("目标可达", { exact: true }).waitFor({ timeout: 60000 });
const downloaded = page.waitForEvent("download");
await page.getByRole("button", { name: "保存场景", exact: true }).click();
await (await downloaded).saveAs(resolve("../artifacts/scene.json"));
await page
  .getByText("打开场景", { exact: true })
  .locator("input")
  .setInputFiles(resolve("../artifacts/scene.json"));
await page.waitForFunction(
  () =>
    document.querySelector('select[aria-label="姿态保持方式"]')?.value ===
    "tracking",
);
await page
  .getByText("导入轨道 JSON / CSV", { exact: true })
  .locator("input")
  .setInputFiles(resolve("../examples/synthetic_orbit.json"));
await page
  .getByText("真实轨道范围内插值，超出覆盖不外推", { exact: false })
  .waitFor({ timeout: 10000 });
await page.getByLabel("北京时间").fill("2027-01-01T00:00");
await page
  .getByRole("alert")
  .filter({ hasText: "超出导入轨道覆盖范围" })
  .waitFor({ timeout: 60000 });
await page.setViewportSize({ width: 390, height: 844 });
await page.waitForFunction(() => document.documentElement.scrollWidth <= 390);
await page.waitForTimeout(400);
await page.screenshot({
  path: resolve("../artifacts/mobile.png"),
  fullPage: true,
});
console.log(
  JSON.stringify({
    pageErrors: errors,
    desktopScreenshot: "digitalQQ/artifacts/desktop.png",
    checks: [
      "WebGL rendering",
      "manual angles",
      "solve",
      "unreachable far target",
      "hold modes",
      "time axis",
      "scene save/load",
      "orbit import",
      "coverage error",
      "mobile layout",
      "PNG download",
      "ZIP report download",
    ],
  }),
);
await browser.close();
if (errors.length) process.exitCode = 1;
