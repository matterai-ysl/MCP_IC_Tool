import asyncio
from types import SimpleNamespace
from concurrent.futures import ThreadPoolExecutor

import pytest

from src.vasp_server import analysis_agent


def test_create_qwen_chat_model_reads_settings(monkeypatch):
    monkeypatch.setattr(analysis_agent.settings, "qwen_api_key", "qwen-key")
    monkeypatch.setattr(
        analysis_agent.settings,
        "qwen_base_url",
        "https://dashscope.aliyuncs.com/compatible-mode/v1",
    )

    captured = {}

    class FakeChatOpenAI:
        def __init__(self, **kwargs):
            captured.update(kwargs)

    def fake_import_module(name):
        assert name == "langchain_openai"
        return SimpleNamespace(ChatOpenAI=FakeChatOpenAI)

    monkeypatch.setattr(analysis_agent.importlib, "import_module", fake_import_module)

    analysis_agent._create_chat_model("qwen3.5-plus")

    assert captured["model"] == "qwen3.5-plus"
    assert captured["api_key"] == "qwen-key"
    assert captured["base_url"] == "https://dashscope.aliyuncs.com/compatible-mode/v1"


def test_create_anthropic_chat_model_reads_settings(monkeypatch):
    monkeypatch.setattr(analysis_agent.settings, "anthropic_api_key", "anthropic-key")

    captured = {}

    class FakeChatAnthropic:
        def __init__(self, **kwargs):
            captured.update(kwargs)

    def fake_import_module(name):
        assert name == "langchain_anthropic"
        return SimpleNamespace(ChatAnthropic=FakeChatAnthropic)

    monkeypatch.setattr(analysis_agent.importlib, "import_module", fake_import_module)

    analysis_agent._create_chat_model("claude-3-7-sonnet")

    assert captured["model"] == "claude-3-7-sonnet"
    assert captured["api_key"] == "anthropic-key"


def test_safe_repl_runs_without_signal_in_worker_thread():
    repl = analysis_agent._SafeREPL({"value": 1})

    with ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(repl.run, "print(value + 1)")
        result = future.result(timeout=5)

    assert "2" in result


def test_run_analysis_passes_recursion_limit(monkeypatch, tmp_path):
    monkeypatch.setattr(analysis_agent.settings, "agent_analysis_max_iterations", 7)
    monkeypatch.setattr(analysis_agent.settings, "agent_analysis_timeout_seconds", 5.0)
    monkeypatch.setattr(analysis_agent, "_make_tools", lambda work_dir: ["tool-a"])
    monkeypatch.setattr(analysis_agent, "_create_chat_model", lambda model_name: object())

    captured = {}

    class FakeAgent:
        async def ainvoke(self, payload, config=None):
            captured["payload"] = payload
            captured["config"] = config
            return {
                "messages": [
                    SimpleNamespace(type="tool", content="listed"),
                    SimpleNamespace(type="ai", content="分析完成"),
                ]
            }

    def fake_import_module(name):
        assert name == "langgraph.prebuilt"
        return SimpleNamespace(create_react_agent=lambda model, tools, prompt: FakeAgent())

    monkeypatch.setattr(analysis_agent.importlib, "import_module", fake_import_module)

    result = asyncio.run(analysis_agent.run_analysis(str(tmp_path), "请总结"))

    assert result["answer"] == "分析完成"
    assert result["steps"] == 1
    assert captured["config"]["recursion_limit"] == 7


def test_run_analysis_times_out(monkeypatch, tmp_path):
    monkeypatch.setattr(analysis_agent.settings, "agent_analysis_max_iterations", 7)
    monkeypatch.setattr(analysis_agent.settings, "agent_analysis_timeout_seconds", 0.01)
    monkeypatch.setattr(analysis_agent, "_make_tools", lambda work_dir: ["tool-a"])
    monkeypatch.setattr(analysis_agent, "_create_chat_model", lambda model_name: object())

    class FakeAgent:
        async def ainvoke(self, payload, config=None):
            await asyncio.sleep(0.05)
            return {"messages": [SimpleNamespace(type="ai", content="分析完成")]}

    def fake_import_module(name):
        assert name == "langgraph.prebuilt"
        return SimpleNamespace(create_react_agent=lambda model, tools, prompt: FakeAgent())

    monkeypatch.setattr(analysis_agent.importlib, "import_module", fake_import_module)

    with pytest.raises(RuntimeError, match="Agent 分析超时"):
        asyncio.run(analysis_agent.run_analysis(str(tmp_path), "请总结"))
