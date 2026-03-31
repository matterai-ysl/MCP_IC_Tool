import types

import pytest

from src.vasp_server import analysis_agent


def test_create_chat_model_uses_qwen_openai_compatible_client(monkeypatch):
    captured = {}

    class FakeChatOpenAI:
        def __init__(self, **kwargs):
            captured.update(kwargs)

    def fake_import_module(name: str):
        if name == "langchain_openai":
            return types.SimpleNamespace(ChatOpenAI=FakeChatOpenAI)
        raise AssertionError(f"unexpected import: {name}")

    monkeypatch.setenv("QWEN_API_KEY", "qwen-key")
    monkeypatch.setenv("QWEN_BASE_URL", "https://qwen.example.com/v1")
    monkeypatch.setattr(analysis_agent.importlib, "import_module", fake_import_module)

    model = analysis_agent._create_chat_model("qwen3.5-plus")

    assert isinstance(model, FakeChatOpenAI)
    assert captured["model"] == "qwen3.5-plus"
    assert captured["api_key"] == "qwen-key"
    assert captured["base_url"] == "https://qwen.example.com/v1"
    assert captured["temperature"] == 0


def test_create_chat_model_uses_anthropic_client_for_claude_models(monkeypatch):
    captured = {}

    class FakeChatAnthropic:
        def __init__(self, **kwargs):
            captured.update(kwargs)

    def fake_import_module(name: str):
        if name == "langchain_anthropic":
            return types.SimpleNamespace(ChatAnthropic=FakeChatAnthropic)
        raise AssertionError(f"unexpected import: {name}")

    monkeypatch.setenv("ANTHROPIC_API_KEY", "anthropic-key")
    monkeypatch.setattr(analysis_agent.importlib, "import_module", fake_import_module)

    model = analysis_agent._create_chat_model("claude-sonnet-4-6")

    assert isinstance(model, FakeChatAnthropic)
    assert captured["model"] == "claude-sonnet-4-6"
    assert captured["api_key"] == "anthropic-key"
    assert captured["temperature"] == 0


def test_create_chat_model_rejects_missing_qwen_env(monkeypatch):
    monkeypatch.delenv("QWEN_API_KEY", raising=False)
    monkeypatch.delenv("QWEN_BASE_URL", raising=False)

    with pytest.raises(ValueError, match="QWEN_API_KEY"):
        analysis_agent._create_chat_model("qwen3.5-plus")
